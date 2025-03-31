import sys
import subprocess
import yaml

REQUIRED_PREFIX = "SYS"
REQUIRED_COVERAGE = 0.5
VALID_TEST_STATUSES = {"implemented", "passed", "ok"}


def load_exported_yaml(prefix: str):
    try:
        result = subprocess.run(
            ["doorstop", "export", prefix, "-y"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            encoding="utf-8"
        )
        return yaml.safe_load(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to export document '{prefix}':\n{e.stderr}", file=sys.stderr)
        sys.exit(1)


def main():
    data = load_exported_yaml(REQUIRED_PREFIX)

    if not isinstance(data, list):
        print(f"ERROR: Exported data for '{REQUIRED_PREFIX}' is not a list.", file=sys.stderr)
        sys.exit(1)

    active_items = [item for item in data if item.get("active") is True]
    if not active_items:
        print(f"ERROR: No active items found in document '{REQUIRED_PREFIX}'.", file=sys.stderr)
        sys.exit(1)

    tested_items = [
        item for item in active_items
        if str(item.get("test", "")).lower() in VALID_TEST_STATUSES
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
