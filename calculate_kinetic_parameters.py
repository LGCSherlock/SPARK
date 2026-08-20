#!/usr/bin/env python3

from math import e
from pathlib import Path
from argparse import ArgumentParser
from scipy.stats import linregress
import pandas as pd
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint


def calc_gamma(path, output, start, end):
    _output = Path(output)
    _output.mkdir(exist_ok=True, parents=True)

    _dir = sorted(Path(path).rglob("CellAnalysis"))
    _barcode = pd.read_csv(_dir[0] / "../GeneMap/Barcode.csv")

    df = pd.DataFrame(
        columns=[
            "gamma",
            "gammar2",
            "lambda",
            "lambdar2",
            "alpha",
            "beta",
            "betar2",
        ]
    )

    print("gamma...")

    # 1. Cytoplasmic RNA translocation parameter gamma
    _all = None
    _index = 0

    for _path in _dir:
        _fn = _path / "CellGene.csv"

        ed = pd.read_csv(_path / "ExpandDistance.csv")
        cg = pd.read_csv(_fn)

        if cg.drop_duplicates().shape != cg.shape:
            print(f"{_fn} has duplicate cell")

        gd = pd.merge(
            ed,
            cg,
            on=("cell", "x", "y"),
        )

        _m = gd.groupby("gene").agg(
            {"dr": ["mean"]}
        )

        _m["file"] = _path.parts[-3]
        _m["time"] = _index

        _index += 2

        _all = pd.concat(
            [_all, _m]
        )

    _all.columns = [
        "dr_mean",
        "file",
        "time",
    ]

    _all.to_csv(
        _output / "Distance.csv"
    )

    plt.figure()
    plt.title("gamma")

    for gene, bc, color, marker in _barcode.values:
        _m = _all.loc[gene]

        _y = _m["dr_mean"]
        _x = _m["time"]

        plt.scatter(
            _x,
            _y,
            c=color,
            marker=marker,
        )

        slope, intercept, r_value, p_value, slope_std_error = linregress(
            _x.iloc[:-1],
            _y.iloc[:-1],
        )

        print(
            f"gene: {gene}, "
            f"slope: {slope}, "
            f"intercept: {intercept}, "
            f"r_value: {r_value}, "
            f"p_value: {p_value}, "
            f"slope_std_error: {slope_std_error}"
        )

        predict_y = (
            intercept
            + slope * _x
        )

        _r2 = (
            0
            if np.isnan(slope)
            else r2_score(
                _y,
                predict_y,
            )
        )

        df.loc[
            gene,
            ["gamma", "gammar2"],
        ] = [
            slope,
            _r2,
        ]

        if not np.isnan(slope):
            plt.plot(
                _x,
                predict_y,
                c=color,
                label=f"{gene}:{_r2:0.3}",
            )

    plt.legend(
        title="Gene",
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
    )

    plt.savefig(
        _output / "Distance.png",
        bbox_inches="tight",
    )

    print("lambda...")

    # 2. Nuclear RNA export parameter lambda
    _all = None
    _index = 0

    for _path in _dir:
        cg = pd.read_csv(
            _path / "CellGene_nuclei.csv"
        )

        ed = pd.read_csv(
            _path / "CellGene.csv"
        )

        ed = ed.groupby("gene").agg(
            {"cell": ["count"]}
        )

        cg = cg.groupby("gene").agg(
            {"cell": ["count"]}
        )

        ed.columns = ["count"]
        cg.columns = ["count_nuclei"]

        gd = pd.merge(
            ed,
            cg,
            on="gene",
        )

        gd["count_cytoplasm"] = (
            gd["count"]
            - gd["count_nuclei"]
        )

        gd["rate_nuclei"] = (
            gd["count_nuclei"]
            / gd["count"]
        )

        gd["file"] = _path.parts[-3]
        gd["time"] = _index

        _index += 2

        _all = pd.concat(
            [_all, gd]
        )

    _all.to_csv(
        _output / "CellGene_Total.csv"
    )

    plt.figure()
    plt.title("lambda")

    for gene, bc, color, marker in _barcode.values:
        _m = _all.loc[gene]

        _y = _m["rate_nuclei"]
        _x = _m["time"]

        plt.scatter(
            _x,
            _y,
            c=color,
            marker=marker,
        )

        slope, intercept, r_value, p_value, slope_std_error = linregress(
            _x.iloc[:-1],
            _y.iloc[:-1],
        )

        print(
            f"gene: {gene}, "
            f"slope: {slope}, "
            f"intercept: {intercept}, "
            f"r_value: {r_value}, "
            f"p_value: {p_value}, "
            f"slope_std_error: {slope_std_error}"
        )

        predict_y = (
            intercept
            + slope * _x
        )

        _r2 = (
            0
            if np.isnan(slope)
            else r2_score(
                _y,
                predict_y,
            )
        )

        df.loc[
            gene,
            ["lambda", "lambdar2"],
        ] = [
            -slope,
            _r2,
        ]

        if not np.isnan(slope):
            plt.plot(
                _x,
                predict_y,
                c=color,
                label=f"{gene}:{_r2:0.3}",
            )

    plt.legend(
        title="Gene",
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
    )

    plt.savefig(
        _output / "CellGene_Total.png",
        bbox_inches="tight",
    )

    print("alpha, beta ...")

    # 3. RNA synthesis parameter alpha and degradation parameter beta
    _all = None
    _index = 0

    # Voxel size used for density normalization
    _ts = 200 * 200 * 350.0

    for _path in _dir:
        cg = pd.read_csv(
            _path / "Matrix.csv"
        )

        _m = pd.DataFrame(
            {
                "sum": cg.sum()[1:].astype(float),
                "count": cg.count()[1:].astype(int),
            }
        )

        _m["desity"] = (
            _m["sum"]
            / (_m["count"] * _ts)
        )

        _m["file"] = _path.parts[-3]
        _m["time"] = _index

        _index += 2

        _all = pd.concat(
            [_all, _m]
        )

    _all.to_csv(
        _output / "GeneMatrix.csv",
        index_label="gene",
    )

    _fig, _axs = plt.subplots(
        len(_barcode),
        1,
        figsize=(20, 10),
        dpi=96,
    )

    _ax_index = 0
    _txy = None

    for gene, bc, color, marker in _barcode.values:
        _m = _all.loc[gene]

        _x = _m["time"]
        _y = _m["desity"]

        # Degradation model:
        # dX(t)/dt = -beta * X(t)
        # X(t) = X0 * exp(-beta * t)

        def model(t, X0, beta):
            return X0 * np.exp(
                -beta * t
            )

        # Degradation model with fixed X0
        def model1(t, beta):
            return X0 * np.exp(
                -beta * t
            )

        # Synthesis model:
        # dX(t)/dt = alpha - beta * X(t)
        # X(t) = (alpha / beta) * (1 - exp(-beta * t))

        def model2(t, alpha, beta):
            return (
                alpha
                / beta
                * (
                    1
                    - np.exp(
                        -beta * t
                    )
                )
            )

        _begin = start - 1

        _end = (
            _x.shape[0] + end
            if end < 0
            else end - 1
            if end <= _x.shape[0]
            else _x.shape[0] - 1
        )

        print(
            f"gene={gene} "
            f"_start={_begin}, "
            f"_end={_end} "
            f"start={start} "
            f"end={end}"
        )

        t_data = _x.iloc[
            _begin:_end + 1
        ]

        X_data = _y.iloc[
            _begin:_end + 1
        ]

        try:
            params, _ = curve_fit(
                model,
                t_data,
                X_data,
            )

        except Exception as e:
            print(
                f"gene={gene} "
                f"error={e} "
                f"continue"
            )
            continue

        X0, beta = params

        fitted_curve_x = np.linspace(
            _x.iloc[0],
            _x.iloc[-1],
            100,
        )

        fitted_curve_y = model(
            fitted_curve_x,
            X0,
            beta,
        )

        print(
            f"{fitted_curve_y[0]} "
            f"?= {model(0, X0, beta)} "
            f"?= {X0} "
            f"!= {_y.iloc[0]}"
        )

        alpha = (
            _y.iloc[0]
            / (
                1
                - np.exp(-beta)
            )
            * beta
        )

        predict_y = model(
            _x,
            X0,
            beta,
        )

        _r2 = r2_score(
            _y.iloc[
                _begin:_end + 1
            ],
            predict_y.iloc[
                _begin:_end + 1
            ],
        )

        df.loc[
            gene,
            [
                "beta",
                "betar2",
                "alpha",
                "X0",
            ],
        ] = [
            beta,
            _r2,
            alpha,
            X0,
        ]

        print(
            f"gene={gene} "
            f"beta={beta} "
            f"alpha={alpha}, "
            f"X0={X0}, "
            f"y[0]={_y.iloc[0]}, "
            f"r2={_r2:0.3}"
        )

        _ax = _axs[
            _ax_index
        ]

        _ax_index += 1

        _ax.set_title(
            f"alpha={alpha} beta={beta}"
        )

        _ax.scatter(
            _x,
            _y,
            c=color,
        )

        # Fitted interval
        fitted_curve_x = np.linspace(
            _x.iloc[_begin],
            _x.iloc[_end],
            100,
        )

        fitted_curve_y = model(
            fitted_curve_x,
            X0,
            beta,
        )

        _txy = pd.concat(
            [
                _txy,
                pd.DataFrame(
                    {
                        "gene": gene,
                        "Tx": fitted_curve_x,
                        "Ty": fitted_curve_y,
                    }
                ),
            ],
            ignore_index=True,
        )

        _ax.plot(
            fitted_curve_x,
            fitted_curve_y,
            color=color,
            label=f"{gene}:{_r2:0.3}",
        )

        _ax.legend(
            title="Gene"
        )

        # Pulse interval
        fitted_curve_x = np.linspace(
            _x.iloc[0],
            _x.iloc[1] // 2,
            30,
        )

        fitted_curve_y = model2(
            fitted_curve_x,
            alpha,
            beta,
        )

        _txy = pd.concat(
            [
                _txy,
                pd.DataFrame(
                    {
                        "gene": gene,
                        "Tx": fitted_curve_x,
                        "Ty": fitted_curve_y,
                    }
                ),
            ],
            ignore_index=True,
        )

        fitted_curve_x -= 1

        _ax.plot(
            fitted_curve_x,
            fitted_curve_y,
            color=color,
        )

        # Extrapolated intervals
        if _begin > 0:
            fitted_curve_x = np.linspace(
                _x.iloc[0],
                _x.iloc[_begin],
                30,
            )

            fitted_curve_y = model(
                fitted_curve_x,
                X0,
                beta,
            )

            _txy = pd.concat(
                [
                    _txy,
                    pd.DataFrame(
                        {
                            "gene": gene,
                            "Tx": fitted_curve_x,
                            "Ty": fitted_curve_y,
                        }
                    ),
                ],
                ignore_index=True,
            )

            _ax.plot(
                fitted_curve_x,
                fitted_curve_y,
                color=color,
                linestyle="--",
            )

        if _end < _x.shape[0] - 1:
            fitted_curve_x = np.linspace(
                _x.iloc[_end],
                _x.iloc[
                    _x.shape[0] - 1
                ],
                30,
            )

            fitted_curve_y = model(
                fitted_curve_x,
                X0,
                beta,
            )

            _txy = pd.concat(
                [
                    _txy,
                    pd.DataFrame(
                        {
                            "gene": gene,
                            "Tx": fitted_curve_x,
                            "Ty": fitted_curve_y,
                        }
                    ),
                ],
                ignore_index=True,
            )

            _ax.plot(
                fitted_curve_x,
                fitted_curve_y,
                color=color,
                linestyle="--",
            )

        fitted_curve_x = np.linspace(
            _x.iloc[1] // 2,
            _x.iloc[1],
            30,
        )

        fitted_curve_y = model2(
            fitted_curve_x,
            alpha,
            beta,
        )

        fitted_curve_x -= 1

        _txy = pd.concat(
            [
                _txy,
                pd.DataFrame(
                    {
                        "gene": gene,
                        "Tx": fitted_curve_x,
                        "Ty": fitted_curve_y,
                    }
                ),
            ],
            ignore_index=True,
        )

        _ax.plot(
            fitted_curve_x,
            fitted_curve_y,
            color=color,
            linestyle="--",
        )

    plt.legend(
        title="Gene",
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
    )

    plt.savefig(
        _output / "GeneMatrix.png",
        bbox_inches="tight",
    )

    df.to_csv(
        _output / "GeneConstant.csv",
        index_label="gene",
    )

    _txy.to_csv(
        _output / "GeneXY.csv"
    )


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Calculate RNA kinetic parameters from CellAnalysis files."
    )

    parser.add_argument(
        "in_path",
        help="Input directory containing CellAnalysis results.",
    )

    parser.add_argument(
        "output_path",
        help="Output directory.",
    )

    parser.add_argument(
        "-s",
        "--start",
        type=int,
        default=2,
        help="Starting sequence index used for fitting (default: %(default)s).",
    )

    parser.add_argument(
        "-e",
        "--end",
        type=int,
        default=-2,
        help="Ending sequence index used for fitting (default: %(default)s).",
    )

    args = parser.parse_args()

    calc_gamma(
        args.in_path,
        args.output_path,
        args.start,
        args.end,
    )
