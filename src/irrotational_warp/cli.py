import argparse
from pathlib import Path

from .viz import plot_slice


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="irrotational_warp")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_plot = sub.add_parser("plot-slice", help="Render a 2D z=0 slice of rho_adm")
    p_plot.add_argument("--rho", type=float, default=10.0, help="Bubble radius (geometric units)")
    p_plot.add_argument("--sigma", type=float, default=5.0, help="Wall sharpness (1/thickness)")
    p_plot.add_argument("--v", type=float, default=1.5, help="Dimensionless v/c")
    p_plot.add_argument("--extent", type=float, default=20.0, help="Half-width of plot domain")
    p_plot.add_argument("--n", type=int, default=301, help="Grid resolution per axis")
    p_plot.add_argument("--out", type=Path, default=Path("results/slice.png"))
    p_plot.add_argument("--json-out", type=Path, default=Path("results/summary.json"))

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "plot-slice":
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        plot_slice(
            rho=args.rho,
            sigma=args.sigma,
            v=args.v,
            extent=args.extent,
            n=args.n,
            out_path=args.out,
            json_out_path=args.json_out,
        )
        return 0

    raise ValueError(f"Unknown command: {args.cmd}")
