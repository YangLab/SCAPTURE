#!/bin/bash

scapture_resolve_path() {
    local caller_dir="${1:-}"
    if [ -n "${SCAPTUREPATH:-}" ] && [ -d "${SCAPTUREPATH%/}/scripts" ]; then
        [[ "$SCAPTUREPATH" != */ ]] && SCAPTUREPATH="${SCAPTUREPATH}/"
        return 0
    fi
    if [ -n "$caller_dir" ] && [ -d "${caller_dir%/}/scripts" ]; then
        SCAPTUREPATH="${caller_dir%/}/"
        return 0
    fi
    local sc_bin sc_real sc_dir
    sc_bin=$(command -v scapture 2>/dev/null || true)
    if [ -n "$sc_bin" ]; then
        sc_real=$(readlink -f "$sc_bin" 2>/dev/null || echo "$sc_bin")
        sc_dir=$(dirname "$sc_real")
        if [ -d "${sc_dir%/}/scripts" ]; then
            SCAPTUREPATH="${sc_dir%/}/"
            return 0
        fi
    fi
    if [ -d "/opt/scapture/scripts" ]; then
        SCAPTUREPATH="/opt/scapture/"
        return 0
    fi
    echo "scapture is not found!" >&2
    return 1
}

scapture_validate_strand() {
    case "${1:-}" in
        forward|reverse|unstranded) return 0 ;;
        *)
            echo "ERROR: Invalid --strand '${1:-}'. Expected one of: forward, reverse, unstranded" >&2
            return 1
            ;;
    esac
}

scapture_featurecounts_strand() {
    scapture_validate_strand "${1:-}" || return 1
    case "${1:-}" in
        reverse) echo 2 ;;
        unstranded) echo 0 ;;
        forward) echo 1 ;;
    esac
}

scapture_samtools_strand_flags() {
    local gene_strand="${1:-}"
    local lib_strand="${2:-}"
    scapture_validate_strand "$lib_strand" || return 1
    if [ "$lib_strand" = "unstranded" ]; then
        echo ""
        return 0
    fi
    if [ "$lib_strand" = "reverse" ]; then
        if [ "$gene_strand" = "+" ]; then
            echo "-f 0x10"
        else
            echo "-F 0x10"
        fi
    else
        if [ "$gene_strand" = "+" ]; then
            echo "-F 0x10"
        else
            echo "-f 0x10"
        fi
    fi
}
