#!/usr/bin/env python3
import argparse, json, gzip, sys
from typing import Dict, Iterable, List, Set, Any, Tuple, Optional, Union

def open_maybe_gzip(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode, encoding=None if "b" in mode else "utf-8")
    return open(path, mode, encoding=None if "b" in mode else "utf-8")

# ---------- dataset parsing (find unlabeled bubble_ids) ----------

def normalize_id(x: Any) -> Optional[int]:
    if isinstance(x, int):
        return x
    if isinstance(x, str):
        try:
            return int(x)
        except ValueError:
            return None
    return None

def read_unlabeled_ids(dataset_jsonl: str, id_field: str = "bubble_id") -> Set[int]:
    """
    Reads your GNN dataset .jsonl and returns the set of bubble_ids that are *unlabeled*.
    A record is considered labeled if any of these are present and not null:
      - label_path
      - label
      - choice_label
    """
    unlabeled: Set[int] = set()
    labeled: Set[int] = set()

    with open_maybe_gzip(dataset_jsonl, "rt") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception as e:
                print(f"[warn] skipping bad JSON line {line_no}: {e}", file=sys.stderr)
                continue

            if id_field not in rec:
                print(f"[warn] line {line_no} missing '{id_field}', skipping", file=sys.stderr)
                continue

            bid = normalize_id(rec[id_field])
            if bid is None:
                print(f"[warn] line {line_no} has non-integer {id_field}={rec[id_field]!r}, skipping", file=sys.stderr)
                continue

            is_labeled = (
                rec.get("label_path", None) is not None
                or rec.get("label", None) is not None
                or rec.get("choice_label", None) is not None
            )

            if is_labeled:
                labeled.add(bid)
            else:
                unlabeled.add(bid)

    return unlabeled - labeled

# ---------- BubbleGun shape helpers ----------

def looks_like_bubble(x: Any) -> bool:
    return isinstance(x, dict) and "id" in x and "ends" in x

def iter_bubbles_any(obj: Any) -> Iterable[Dict[str, Any]]:
    """
    Iterate every bubble dict anywhere we can find a 'bubbles' list.
    """
    if isinstance(obj, dict):
        if "bubbles" in obj and isinstance(obj["bubbles"], list):
            for b in obj["bubbles"]:
                if looks_like_bubble(b):
                    yield b
        else:
            for v in obj.values():
                yield from iter_bubbles_any(v)
    elif isinstance(obj, list):
        for v in obj:
            if looks_like_bubble(v):
                yield v
            else:
                yield from iter_bubbles_any(v)

def filter_bubbles_list(lst: List[Dict[str, Any]], wanted_ids: Set[int], seen: Set[int]) -> List[Dict[str, Any]]:
    out = []
    for b in lst:
        if not looks_like_bubble(b):
            continue
        bid = normalize_id(b.get("id"))
        if bid is None:
            continue
        if bid in wanted_ids:
            out.append(b)
            seen.add(bid)
    return out

def mirror_shape_and_filter(obj: Any, wanted_ids: Set[int], drop_empty_chains: bool, seen: Set[int]) -> Any:
    """
    Return a new object with the *same top-level shape* as obj, but with bubbles filtered
    to those whose id is in wanted_ids. Preserves other metadata.
    """
    # Case 1: top-level dict containing "bubbles": [...]
    if isinstance(obj, dict) and "bubbles" in obj and isinstance(obj["bubbles"], list):
        new_obj = dict(obj)  # shallow copy
        new_obj["bubbles"] = filter_bubbles_list(obj["bubbles"], wanted_ids, seen)
        return new_obj

    # Case 2: dict-of-chains: {"1": {"bubbles":[...]}, "2": {...}, ...}
    if isinstance(obj, dict):
        out_map: Dict[str, Any] = {}
        for k, v in obj.items():
            if isinstance(v, dict) and "bubbles" in v and isinstance(v["bubbles"], list):
                nv = dict(v)
                nv["bubbles"] = filter_bubbles_list(v["bubbles"], wanted_ids, seen)
                if nv["bubbles"] or not drop_empty_chains:
                    out_map[k] = nv
            else:
                # recurse in case chains are nested oddly
                nv = mirror_shape_and_filter(v, wanted_ids, drop_empty_chains, seen)
                # Keep non-bubble entries as-is (e.g., metadata blocks), unless you want to drop:
                out_map[k] = nv
        return out_map

    # Case 3: list-of-chains or list-of-bubbles
    if isinstance(obj, list) and obj:
        # List of chains (each element has "bubbles")
        if isinstance(obj[0], dict) and "bubbles" in obj[0]:
            out_list = []
            for chain in obj:
                if isinstance(chain, dict) and "bubbles" in chain and isinstance(chain["bubbles"], list):
                    nc = dict(chain)
                    nc["bubbles"] = filter_bubbles_list(chain["bubbles"], wanted_ids, seen)
                    if nc["bubbles"] or not drop_empty_chains:
                        out_list.append(nc)
                else:
                    out_list.append(chain)
            return out_list
        # Flat list of bubbles
        elif looks_like_bubble(obj[0]):
            return filter_bubbles_list(obj, wanted_ids, seen)

    # Fallback: return unchanged (no bubbles found here)
    return obj

# ---------- glue ----------

def main():
    ap = argparse.ArgumentParser(description="Extract unlabeled bubbles and emit JSON mirroring BubbleGun's original container shape.")
    ap.add_argument("--dataset", required=True, help="Path to your GNN dataset .jsonl (has bubble_id, label_path/choice_label/label)")
    ap.add_argument("--bubblegun", required=True, help="Path to BubbleGun JSON (chains+bubbles)")
    ap.add_argument("--out", required=True, help="Output JSON path (will preserve the input's top-level shape)")
    ap.add_argument("--id-field", default="bubble_id", help="Dataset id field name (default: bubble_id)")
    ap.add_argument("--keep-empty-chains", action="store_true", help="Keep chains that end up with zero bubbles after filtering")
    args = ap.parse_args()

    wanted_ids = read_unlabeled_ids(args.dataset, id_field=args.id_field)
    print(f"[info] dataset unlabeled ids: {len(wanted_ids)}", file=sys.stderr)

    with open_maybe_gzip(args.bubblegun, "rt") as f:
        src = json.load(f)

    seen: Set[int] = set()
    filtered = mirror_shape_and_filter(src, wanted_ids, drop_empty_chains=(not args.keep_empty_chains), seen=seen)

    missing = wanted_ids - seen
    if missing:
        mlist = sorted(list(missing))[:20]
        print(f"[note] {len(missing)} requested ids not present in BubbleGun JSON (first few: {mlist})", file=sys.stderr)

    with open_maybe_gzip(args.out, "wt") as w:
        json.dump(filtered, w, indent=2)

    # quick counts for sanity (how many bubbles are in the *output*)
    out_bubbles = sum(1 for _ in iter_bubbles_any(filtered))
    print(f"[done] wrote {out_bubbles} bubbles to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
