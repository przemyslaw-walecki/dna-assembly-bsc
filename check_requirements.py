import sys
import yaml

REQUIRED_PREFIX = "SYS"
REQUIRED_COVERAGE = 0.5
VALID_TEST_STATUSES = {"implemented", "passed", "ok"}
EXPORT_FILE = "requirements.yaml"

def main():
    try:
        with open(EXPORT_FILE, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"ERROR: Exported file '{EXPORT_FILE}' not found.", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"ERROR: Failed to parse YAML: {e}", file=sys.stderr)
        sys.exit(1)

    if REQUIRED_PREFIX not in data:
        print(f"ERROR: No document with prefix '{REQUIRED_PREFIX}' in exported data.", file=sys.stderr)
        sys.exit(1)

    document = data[REQUIRED_PREFIX]

    active_items = [item for item in document.values() if item.get("active", False)]
    if not active_items:
        print(f"ERROR: No active items in document '{REQUIRED_PREFIX}'.", file=sys.stderr)
        sys.exit(1)

    tested_items = [
        item for item in active_items
        if isinstance(item.get("test", ""), str) and item["test"].lower() in VALID_TEST_STATUSES
    ]

    total = len(active_items)
    tested = len(tested_items)
    coverage = tested / total

    print(f"{REQUIRED_PREFIX} coverage: {tested}/{total} ({coverage:.1%})")

    if coverage < REQUIRED_COVERAGE:
        print(f"ERROR: Test coverage below threshold ({REQUIRED_COVERAGE:.0%}).", file=sys.stderr)
        sys.exit(1)

    print("SUCCESS: Requirement test coverage OK.")
    sys.exit(0)

if __name__ == "__main__":
    main()
