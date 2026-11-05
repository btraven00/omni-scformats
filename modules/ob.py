"""Parse the omnibenchmark 0.8 module CLI:
    --output id=PATH (repeated)  --<input_id> PATH  --<param> VALUE  --output_dir/--name
Gather inputs repeat their flag, so every value is kept as a list."""
import sys


def parse(argv=None):
    args, outputs, key = {}, {}, None
    for tok in argv if argv is not None else sys.argv[1:]:
        if tok.startswith("--"):
            key = tok[2:]
        elif key == "output":
            oid, path = tok.split("=", 1)
            outputs[oid] = path
        elif key is not None:
            args.setdefault(key, []).append(tok)
    one = {k: v[0] for k, v in args.items()}
    return one, args, outputs


if __name__ == "__main__":
    one, many, out = parse(["--output", "a.b=/x=y", "--rds", "/r", "--perf", "1", "--perf", "2"])
    assert out == {"a.b": "/x=y"} and one["rds"] == "/r" and many["perf"] == ["1", "2"], (one, many, out)
    print("ok")
