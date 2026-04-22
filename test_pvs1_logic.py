#!/usr/bin/env python3
"""
Standalone test for PVS1_general NULL variant logic.
Tests the core decision tree without requiring module imports.
"""

from enum import Enum


class VARTYPE(Enum):
    STOP_GAINED = "stop_gained"
    FRAMESHIFT_VARIANT = "frameshift_variant"
    START_LOST = "start_lost"
    SPLICE_DONOR = "splice_donor"
    SPLICE_ACCEPTOR = "splice_acceptor"
    MISSENSE_VARIANT = "missense_variant"
    SYNONYMOUS_VARIANT = "synonymous_variant"


# NULL_LOF set from var_type.py
NULL_LOF = {
    VARTYPE.STOP_GAINED,
    VARTYPE.FRAMESHIFT_VARIANT,
    VARTYPE.SPLICE_DONOR,
    VARTYPE.SPLICE_ACCEPTOR,
    VARTYPE.START_LOST,
}


def is_null_variant(var_types):
    """Check if variant is NULL/LOF type."""
    for var_type in var_types:
        if var_type in NULL_LOF:
            return True
    return False


def get_evidence_strength_frameshift(is_NMD, is_critical_region, diff_len_pct, threshold=10.0):
    """
    Calculate PVS1 evidence strength for frameshift/null variants.

    Logic:
    - NMD → Very Strong
    - Not NMD + Critical region → Strong
    - Not NMD + Not critical region → Moderate
    """
    if is_NMD:
        return "VERY_STRONG"

    if is_critical_region:
        return "STRONG"

    return "MODERATE"


def get_evidence_strength_start_loss(has_alt_start, is_critical_region):
    """
    Calculate PVS1 evidence strength for start loss variants.

    Logic (original ACMG):
    - No alt start codon → Moderate
    - Alt start + Critical region → Moderate
    - Alt start + Non-critical region → Supporting
    """
    if not has_alt_start:
        return "MODERATE"

    if is_critical_region:
        return "MODERATE"

    return "SUPPORTING"


def test_null_variant_detection():
    """Test NULL variant type detection."""
    print("\n=== Test 1: NULL Variant Detection ===")

    test_cases = [
        ([VARTYPE.STOP_GAINED], True, "stop_gained"),
        ([VARTYPE.FRAMESHIFT_VARIANT], True, "frameshift"),
        ([VARTYPE.START_LOST], True, "start_lost"),
        ([VARTYPE.SPLICE_DONOR], True, "splice_donor"),
        ([VARTYPE.SPLICE_ACCEPTOR], True, "splice_acceptor"),
        ([VARTYPE.MISSENSE_VARIANT], False, "missense"),
        ([VARTYPE.SYNONYMOUS_VARIANT], False, "synonymous"),
    ]

    for var_types, expected_null, name in test_cases:
        result = is_null_variant(var_types)
        status = "✓" if result == expected_null else "✗"
        print(f"  {status} {name}: is_null={result} (expected={expected_null})")
        assert result == expected_null, f"Failed for {name}"

    print("  ✓ All NULL variant detection tests passed!")


def test_frameshift_evidence_strength():
    """Test frameshift PVS1 evidence strength."""
    print("\n=== Test 2: Frameshift PVS1 Evidence Strength ===")

    test_cases = [
        # (is_NMD, is_critical, diff_len_pct, expected_strength)
        (True, False, 50.0, "VERY_STRONG"),
        (True, True, 50.0, "VERY_STRONG"),
        (False, True, 50.0, "STRONG"),
        (False, False, 50.0, "MODERATE"),
        (False, False, 5.0, "MODERATE"),
    ]

    for is_NMD, is_critical, diff_len, expected in test_cases:
        result = get_evidence_strength_frameshift(is_NMD, is_critical, diff_len)
        status = "✓" if result == expected else "✗"
        nmd_str = "NMD" if is_NMD else "NoNMD"
        crit_str = "Crit" if is_critical else "NonCrit"
        print(f"  {status} {nmd_str}+{crit_str}+{diff_len}% → {result} (expected={expected})")
        assert result == expected, f"Failed: {nmd_str}+{crit_str} → {result}"

    print("  ✓ All frameshift evidence strength tests passed!")


def test_start_loss_evidence_strength():
    """Test start loss PVS1 evidence strength."""
    print("\n=== Test 3: Start Loss PVS1 Evidence Strength ===")

    test_cases = [
        # (has_alt_start, is_critical, expected_strength)
        (False, False, "MODERATE"),
        (False, True, "MODERATE"),
        (True, True, "MODERATE"),
        (True, False, "SUPPORTING"),
    ]

    for has_alt, is_crit, expected in test_cases:
        result = get_evidence_strength_start_loss(has_alt, is_crit)
        status = "✓" if result == expected else "✗"
        alt_str = "Alt" if has_alt else "NoAlt"
        crit_str = "Crit" if is_crit else "NonCrit"
        print(f"  {status} {alt_str}+{crit_str} → {result} (expected={expected})")
        assert result == expected, f"Failed: {alt_str}+{crit_str} → {result}"

    print("  ✓ All start loss evidence strength tests passed!")


def show_summary():
    """Show PVS1 decision summary."""
    print("\n" + "=" * 60)
    print("PVS1_general Decision Summary")
    print("=" * 60)

    print("\n[Frameshift/Stop_gained/Splice]")
    print("  NULL variant + NMD                    → VERY_STRONG")
    print("  NULL variant + NOT NMD + Critical     → STRONG")
    print("  NULL variant + NOT NMD + Non-critical → MODERATE")

    print("\n[Start Loss]")
    print("  No alt start codon                    → MODERATE")
    print("  Alt start + Critical region           → MODERATE")
    print("  Alt start + Non-critical region       → SUPPORTING")

    print("\n[Splice Variant]")
    print("  Splice alteration + NMD               → VERY_STRONG")
    print("  Splice alteration + Critical         → STRONG")
    print("  Splice alteration + Non-critical     → MODERATE")
    print("  No splice alteration predicted       → NOT APPLICABLE")


def run_tests():
    """Run all tests."""
    print("=" * 60)
    print("PVS1_general Logic Tests")
    print("=" * 60)

    try:
        test_null_variant_detection()
        test_frameshift_evidence_strength()
        test_start_loss_evidence_strength()
        show_summary()

        print("\n" + "=" * 60)
        print("ALL TESTS PASSED!")
        print("=" * 60)
        return True

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        return False


if __name__ == "__main__":
    success = run_tests()
    exit(0 if success else 1)
