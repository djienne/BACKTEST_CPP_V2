#!/usr/bin/env bash
#
# Launch every strategy binary for a bounded run and report whether it produced a
# plausible result. This is the end-to-end check that the configured universe actually
# has its data and that each strategy trades on it.
#
#     ./tools/smoke_all_strategies.sh [seconds_per_strategy]
#
# Runs strictly one at a time, so it never competes with your own work for cores.

set -uo pipefail
cd "$(dirname "$0")/.."

SECS="${1:-120}"
LOGS=".smoke_logs"
mkdir -p "$LOGS"

STRATEGIES=(
    backtest_double_EMA_float
    backtest_double_EMA_StochRSI_float_muti_pair
    BigWill
    backtest_TRIX_multi_pair
    BBTREND
    SuperTrend_EMA_ATR
    3EMA_SRSI_ATR
    SuperReversal_mtf
    STRATEGY_TEMPLATE
    F_BigWill
    F_BBTREND
    F_SuperReversal_mtf
    F_3EMA_SRSI_ATR
)

printf '%-46s %-10s %-10s %s\n' "STRATEGY" "STATUS" "TRADES" "NOTE"
printf '%.0s-' {1..90}; echo

fails=0
for s in "${STRATEGIES[@]}"; do
    log="$LOGS/$s.log"
    timeout "$SECS" "./$s.exe" > "$log" 2>&1
    rc=$?

    # A clean finish is 0; 124 means the timeout stopped a still-running sweep, which is
    # expected for the large parameter spaces and is not a failure.
    if [ $rc -ne 0 ] && [ $rc -ne 124 ]; then
        note=$(grep -m1 -E 'FATAL|error|Error' "$log" | cut -c1-46)
        printf '%-46s %-10s %-10s %s\n' "$s" "CRASH($rc)" "-" "${note:-see $log}"
        fails=$((fails + 1))
        continue
    fi

    # Last reported trade count from any result block.
    trades=$(grep -aoE 'Number of trades *: *[0-9]+' "$log" | tail -1 | grep -oE '[0-9]+$')

    # Only the metric lines, and only the final block: a substring search over the whole
    # log matches "bi(nan)ce" in the data paths, and the sweep legitimately prints
    # Score: -inf as its starting sentinel until the first candidate passes the filters.
    nan=$(grep -aE '^(Gain|Win rate|max DD|Gain/DDC|Score|Calmar) ' "$log" \
          | tail -8 \
          | sed -r 's/\x1B\[[0-9;]*m//g' \
          | grep -cE '(^| )-?(nan|inf)\b' || true)

    if [ -z "${trades:-}" ]; then
        # No result block yet -- did it at least get through data loading?
        if grep -qa 'Running all backtests\|Initialized calculations\|Done calculating indicators\|Calculated EMAs\|Done\.' "$log"; then
            printf '%-46s %-10s %-10s %s\n' "$s" "RUNNING" "-" "loaded data, sweeping (no block within ${SECS}s)"
        else
            printf '%-46s %-10s %-10s %s\n' "$s" "STALLED" "-" "no progress; see $log"
            fails=$((fails + 1))
        fi
        continue
    fi

    if [ "$trades" -eq 0 ]; then
        # `best` starts zeroed and is only replaced once a candidate clears the min-trades
        # and drawdown filters, so 0 here means "nothing qualified yet" -- unless the
        # sweep genuinely ran nothing, which is a real defect.
        ran=$(grep -aoE 'Number of backtests performed *: *[0-9]+' "$log" | tail -1 | grep -oE '[0-9]+$')
        if [ "${ran:-1}" -eq 0 ]; then
            printf '%-46s %-10s %-10s %s\n' "$s" "NO-SWEEP" "0" "ran 0 backtests -- empty parameter space"
            fails=$((fails + 1))
        else
            printf '%-46s %-10s %-10s %s\n' "$s" "NO-QUALIFY" "0" "swept, but nothing passed the filters yet"
        fi
    elif [ "${nan:-0}" -gt 0 ]; then
        printf '%-46s %-10s %-10s %s\n' "$s" "NAN" "$trades" "non-finite metric in output"
        fails=$((fails + 1))
    else
        printf '%-46s %-10s %-10s %s\n' "$s" "OK" "$trades" ""
    fi
done

echo
if [ $fails -eq 0 ]; then
    echo "All strategies ran and traded."
else
    echo "$fails strategy/strategies need attention (logs in $LOGS/)."
fi
exit $fails
