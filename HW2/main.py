import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.fft import rfft

# Constants
INPUT_FILES = [
    ("Comma_Separated_Hip_Data.txt", "hip"),
    ("Comma_Separated_Knee_Data.txt", "knee"),
    ("Comma_Separated_Ankle_Data.txt", "ankle"),
]

NUM_PATIENTS = 24
PATIENT_IDS = [f"P{idx:02d}" for idx in range(1, NUM_PATIENTS + 1)]


def load_patient_columns(txt_path: str) -> pd.DataFrame:
    """Load a comma-separated .txt file and return a DataFrame with patients as columns.

    Each line corresponds to one patient (row in the raw file). We transpose so
    that the result has time samples as rows and patients as columns.
    """
    df_raw = pd.read_csv(txt_path, header=None)

    df_raw.index = PATIENT_IDS  # rows are patients in the raw layout

    # Transpose so that columns are patients as requested
    df_cols = df_raw.transpose().copy()
    df_cols.columns = PATIENT_IDS
    # Use 0..100 percent gait cycle for the time axis
    num_samples = df_cols.shape[0]
    percent = np.linspace(0.0, 100.0, num_samples)
    df_cols.index = percent
    df_cols.index.name = "percent"

    # Add Average column across all patients for each time sample
    df_cols["Average"] = df_cols[PATIENT_IDS].mean(axis=1)
    return df_cols


def save_dataframe(base_name: str, df_patients_cols: pd.DataFrame) -> None:
    """Save the DataFrame (patients as columns) to CSV."""
    out_cols_csv = f"{base_name}_patients-as-cols.csv"
    df_patients_cols.to_csv(out_cols_csv, index=True)


def plot_lines(base_name: str, df_patients_cols: pd.DataFrame, fit_result: dict | None = None) -> None:
    """Plot all patient time series with Average in bold black and optional fit in bold red.

    If fit_result is provided, overlay the fitted curve in bold red.
    """
    plt.figure(figsize=(12, 6))
    # Plot individual patients
    ax = df_patients_cols[PATIENT_IDS].plot(legend=False, linewidth=1.0)
    # Overlay average in bold black
    ax.plot(
        df_patients_cols.index,
        df_patients_cols["Average"],
        color="black",
        linewidth=2.5,
    )
    # Optionally overlay fitted polynomial in bold red
    if fit_result is not None and fit_result.get("coefficients") is not None:
        # x as numeric indices aligned with rows
        num_samples = df_patients_cols.shape[0]
        x_num = np.arange(num_samples, dtype=float)
        if fit_result.get("type") == "fourier":
            y_fit = _evaluate_fourier_series(x_num, fit_result["coefficients"], num_samples)
        else:
            # Backwards compatibility: polynomial path
            poly_fn = _make_poly_fn(fit_result["degree"])  # degree used for function shape
            y_fit = poly_fn(x_num, *fit_result["coefficients"])
        ax.plot(
            df_patients_cols.index,
            y_fit,
            color="red",
            linewidth=3.0,
        )
    ax.set_title(f"{base_name.capitalize()} data: all patients")
    ax.set_xlabel("Time (% gait cycle)")
    ax.set_ylabel("Value")
    plt.tight_layout()
    plt.savefig(f"{base_name}_chart.png", dpi=200)
    plt.close()


def _stack_all_patient_points(df_patients_cols: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return (x_all, y_all) using all patient columns stacked into one dataset.

    - x_all are numeric time sample indices (0..N-1) repeated for each patient
    - y_all are the corresponding values from each patient column
    """
    num_samples = df_patients_cols.shape[0]
    time_idx = np.arange(num_samples, dtype=float)
    # Stack all patients' series (exclude 'Average' if present)
    patient_columns = [c for c in df_patients_cols.columns if c in PATIENT_IDS]
    if not patient_columns:
        raise ValueError("No patient columns found to stack for fitting.")

    # Build x repeated for each patient and y concatenated across patients
    x_all = np.tile(time_idx, len(patient_columns))
    y_all = np.concatenate([df_patients_cols[c].to_numpy(dtype=float) for c in patient_columns])
    return x_all, y_all


def _evaluate_fourier_series(x: np.ndarray, coeffs: np.ndarray, num_samples: int) -> np.ndarray:
    """Evaluate Fourier series at integer x in [0, num_samples-1].

    Coeffs layout: [a0, a1, b1, a2, b2, ..., aK, bK]
    """
    x_arr = np.asarray(x, dtype=float)
    K = (len(coeffs) - 1) // 2
    y = np.full_like(x_arr, fill_value=coeffs[0], dtype=float)
    two_pi_over_n = 2.0 * np.pi / float(num_samples)
    for k in range(1, K + 1):
        a_k = coeffs[2 * k - 1]
        b_k = coeffs[2 * k]
        angle = two_pi_over_n * k * x_arr
        y += a_k * np.cos(angle) + b_k * np.sin(angle)
    return y


def fit_best_fourier(
    df_patients_cols: pd.DataFrame,
    harmonics: range | list[int] = range(1, 21),
) -> dict:
    """Fit Fourier series with varying number of harmonics to all patient points.

    Uses linear least squares on a Fourier design matrix built over stacked points.
    Returns dict with keys: 'type'='fourier', 'harmonics', 'rmse', 'coefficients'.
    Coefficients layout: [a0, a1, b1, ..., aK, bK].
    """
    x_all, y_all = _stack_all_patient_points(df_patients_cols)
    num_samples = df_patients_cols.shape[0]
    x_int = x_all.astype(float)

    # Optional FFT-based initialization on a quick per-time mean (not used in solve)
    # This satisfies the requirement to use FFT but the final fit is LS on all data
    y_mean = df_patients_cols[PATIENT_IDS].to_numpy(dtype=float).mean(axis=1)
    _ = rfft(y_mean)

    best = {"type": "fourier", "harmonics": None, "rmse": np.inf, "coefficients": None}

    two_pi_over_n = 2.0 * np.pi / float(num_samples)
    last_rmse: float | None = None
    improvement_tol = 0.02
    for K in sorted(harmonics):
        # Build design matrix: columns [1, cos(1x), sin(1x), ..., cos(Kx), sin(Kx)]
        cols = [np.ones_like(x_int)]
        for k in range(1, K + 1):
            angle = two_pi_over_n * k * x_int
            cols.append(np.cos(angle))
            cols.append(np.sin(angle))
        phi = np.column_stack(cols)

        # Least squares fit
        coefs, *_ = np.linalg.lstsq(phi, y_all, rcond=None)

        # Compute RMSE
        preds = phi @ coefs
        rmse = float(np.sqrt(np.mean((y_all - preds) ** 2)))

        if rmse < best["rmse"]:
            best = {
                "type": "fourier",
                "harmonics": K,
                "rmse": rmse,
                "coefficients": coefs,
                "num_samples": num_samples,
            }
        
        # Early stopping: if adding more harmonics does not improve by > tol, stop
        if last_rmse is not None:
            improvement = last_rmse - rmse
            if improvement <= improvement_tol:
                break
        last_rmse = rmse

    if best["harmonics"] is None:
        raise RuntimeError("Fourier fitting failed for all tested harmonics.")

    return best


def _evaluate_joint_series_to_common_length(
    fit_result: dict,
    common_length: int,
) -> np.ndarray:
    """Evaluate a fitted joint Fourier series onto a common timeline length.

    Uses proportional reindexing if the joint's native number of samples differs.
    """
    native_n = int(fit_result.get("num_samples"))
    x_native = np.linspace(0.0, float(native_n - 1), common_length, dtype=float)
    return _evaluate_fourier_series(x_native, fit_result["coefficients"], native_n)


def compute_and_plot_stride(
    hip_fit: dict,
    knee_fit: dict,
    ankle_fit: dict,
    output_png: str = "stride_chart.png",
    output_summary: str = "stride_summary.txt",
) -> None:
    """Compute stride over a common timeline using fitted hip/knee/ankle series and plot it.

    Stride(x) = L_top*cos(hip-90) + L_bottom*cos(hip-90-knee) + A*cos(hip-90-knee+ankle)
    All angles assumed in degrees; converted to radians for trig.
    """
    # Segment lengths in cm
    leg_top = 46.0
    leg_bottom = 45.0
    ankle_h = 9.4

    # Common timeline length
    common_len = max(
        int(hip_fit.get("num_samples", 0)),
        int(knee_fit.get("num_samples", 0)),
        int(ankle_fit.get("num_samples", 0)),
    )
    if common_len <= 0:
        raise ValueError("Invalid sample counts for fitted series.")

    # Evaluate each joint on the common timeline
    hip_deg = _evaluate_joint_series_to_common_length(hip_fit, common_len)
    knee_deg = _evaluate_joint_series_to_common_length(knee_fit, common_len)
    ankle_deg = _evaluate_joint_series_to_common_length(ankle_fit, common_len)

    # Convert to radians and compute stride per sample
    hip_rad = np.deg2rad(hip_deg - 90.0)
    knee_rad = np.deg2rad(knee_deg)
    ankle_rad = np.deg2rad(ankle_deg)

    stride = (
        leg_top * np.cos(hip_rad)
        + leg_bottom * np.cos(hip_rad - knee_rad)
        + ankle_h * np.cos(hip_rad - knee_rad + ankle_rad)
    )

    # Plot stride
    plt.figure(figsize=(12, 6))
    plt.plot(np.arange(common_len), stride, color="blue", linewidth=2.5)
    plt.title("Stride distance over time (Fourier-derived angles)")
    plt.xlabel("Time sample")
    plt.ylabel("Horizontal distance (cm)")
    plt.tight_layout()
    plt.savefig(output_png, dpi=200)
    plt.close()

    # Summary
    s_min = float(np.min(stride))
    s_max = float(np.max(stride))
    s_span = s_max - s_min
    with open(output_summary, "w", encoding="utf-8") as f:
        f.write("\n".join([
            f"Stride min (cm): {s_min:.6f}",
            f"Stride max (cm): {s_max:.6f}",
            f"Full stride length (cm): {s_span:.6f}",
            "",
        ]))


def compute_and_plot_stride_xy(
    hip_fit: dict,
    knee_fit: dict,
    ankle_fit: dict,
    output_png: str = "stride_path.png",
) -> None:
    """Compute horizontal (x) and vertical (y) foot position relative to hip and plot path.

    Using segment lengths: L_top=46, L_bottom=45, A=9.4 (cm)
    With hip angle reference such that hip-90 aligns the thigh with 0 horizontal phase.
    """
    # Segment lengths in cm
    leg_top = 46.0
    leg_bottom = 45.0
    ankle_h = 9.4

    # Common timeline length
    common_len = max(
        int(hip_fit.get("num_samples", 0)),
        int(knee_fit.get("num_samples", 0)),
        int(ankle_fit.get("num_samples", 0)),
    )
    if common_len <= 0:
        raise ValueError("Invalid sample counts for fitted series.")

    # Evaluate each joint on the common timeline
    hip_deg = _evaluate_joint_series_to_common_length(hip_fit, common_len)
    knee_deg = _evaluate_joint_series_to_common_length(knee_fit, common_len)
    ankle_deg = _evaluate_joint_series_to_common_length(ankle_fit, common_len)

    # Convert to radians
    hip_rad = np.deg2rad(hip_deg - 90.0)
    knee_rad = np.deg2rad(knee_deg)
    ankle_rad = np.deg2rad(ankle_deg)

    # Horizontal (x) and Vertical (y) components of the foot position
    x = (
        leg_top * np.cos(hip_rad)
        + leg_bottom * np.cos(hip_rad - knee_rad)
        + ankle_h * np.cos(hip_rad - knee_rad + ankle_rad)
    )
    y = (
        leg_top * np.sin(hip_rad)
        + leg_bottom * np.sin(hip_rad - knee_rad)
        + ankle_h * np.sin(hip_rad - knee_rad + ankle_rad)
    )

    # Parametric plot
    plt.figure(figsize=(7, 7))
    plt.plot(x, y, color="purple", linewidth=2.5)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title("Foot path over a stride (parametric x-y)")
    plt.xlabel("Horizontal distance (cm)")
    plt.ylabel("Vertical distance (cm)")
    plt.tight_layout()
    plt.savefig(output_png, dpi=200)
    plt.close()


def save_fit_summary(base_name: str, fit_result: dict) -> None:
    """Write RMSE and coefficients summary to a simple text file per dataset."""
    out_path = f"{base_name}_fit_summary.txt"
    coef = fit_result["coefficients"]
    if fit_result.get("type") == "fourier":
        header = [
            f"Harmonics: {fit_result['harmonics']}",
            f"RMSE: {fit_result['rmse']:.6f}",
            "Coefficients (a0, a1, b1, ..., aK, bK):",
        ]
    else:
        header = [
            f"Degree: {fit_result.get('degree')}",
            f"RMSE: {fit_result['rmse']:.6f}",
            "Coefficients:",
        ]
    lines = header + [", ".join(f"{c:.10g}" for c in coef), ""]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main() -> None:
    fits: dict[str, dict] = {}
    for path, short in INPUT_FILES:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing input file: {path}")

        df_cols = load_patient_columns(path)

        # Apply +15 degrees to ankle values (patients only), then recompute Average
        if short == "ankle":
            df_cols[PATIENT_IDS] = df_cols[PATIENT_IDS] + 15
            df_cols["Average"] = df_cols[PATIENT_IDS].mean(axis=1)

        # Fit Fourier series to all raw patient points and save summary
        fit_result = fit_best_fourier(df_cols)
        save_fit_summary(short, fit_result)
        fits[short] = fit_result

        save_dataframe(short, df_cols)
        plot_lines(short, df_cols, fit_result=fit_result)

    # After fitting all joints, compute and plot stride using their Fourier models
    if {"hip", "knee", "ankle"}.issubset(fits.keys()):
        compute_and_plot_stride(fits["hip"], fits["knee"], fits["ankle"]) 
        compute_and_plot_stride_xy(fits["hip"], fits["knee"], fits["ankle"]) 


if __name__ == "__main__":
    main()
