### ME139 HW2 — Data Import and Plotting (macOS/Linux)

This project reads three comma-separated `.txt` files (`Comma_Separated_Hip_Data.txt`, `Comma_Separated_Knee_Data.txt`, `Comma_Separated_Ankle_Data.txt`) into pandas DataFrames where each column is a patient (`P01`..`P24`), computes an `Average` column per time sample, and saves both CSVs and charts.

### Prerequisites
- macOS or Linux with bash/zsh
- Python 3.9+ (recommended `python3` on your PATH)

### 1) Set up a virtual environment (bash/zsh)
```bash
cd "C:\\Users\\Luke\\Documents\\Classes\\ME139\\HW2"  # If on macOS/Linux, cd to your project folder path

# Create a fresh venv in the project (named .venv)
python3 -m venv .venv

# Activate the venv
source .venv/bin/activate

# Upgrade pip and install required packages
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Note: The `venv/` in this repo is Windows-specific and won’t be usable on macOS/Linux. Use `.venv` as shown above on macOS/Linux.

### 2) Run the script
```bash
python main.py
```

### 3) Outputs
For each input file you’ll get:
- CSV with patients as columns and an `Average` column:
  - `hip_patients-as-cols.csv`
  - `knee_patients-as-cols.csv`
  - `ankle_patients-as-cols.csv` (ankle values include a +15 offset, and `Average` reflects that)
- PNG chart with all patients plotted, the `Average` line in bold black, and the red fitted curve:
  - `hip_chart.png`, `knee_chart.png`, `ankle_chart.png`
- Additional stride plots and summary:
  - `stride_chart.png`, `stride_path.png`, `stride_summary.txt`

### 3b) Merge the three charts into one PDF (Pandoc)

Ensure Pandoc is installed:
- macOS (Homebrew):
  ```bash
  brew install pandoc
  ```
- Windows: download the installer from `https://pandoc.org/installing.html` and add Pandoc to PATH if needed.

From the project directory, run:
```bash
# macOS/Linux (bash/zsh)
pandoc hip_chart.png knee_chart.png ankle_chart.png -V geometry:margin=0.5in -o charts.pdf

# Windows PowerShell / cmd
pandoc hip_chart.png knee_chart.png ankle_chart.png -V geometry:margin=0.5in -o charts.pdf
```

Notes:
- The order of images determines page order in `charts.pdf`.
- If Pandoc prompts for a LaTeX engine, install a minimal TeX distribution (e.g., TinyTeX) or use `-V papersize:letter` as needed.

### 4) Deactivate the virtual environment
```bash
deactivate
```

### Windows setup notes
- Use the included `venv\` only on Windows. On macOS/Linux, create `.venv` as shown.
- On Windows PowerShell you can install via:
```powershell
python -m venv .venv
. .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Optional: Re-create requirements file
```bash
pip freeze > requirements.txt
```
Then others can do:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


