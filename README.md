# CTFFIND3 Parallel Processing & CTF Fit Resolution

Shell scripts to run **CTFFIND3** in parallel from the terminal and to calculate **CTF fit resolution** from CTFFIND3 power spectrum and log files.

---

## Requirements

- Linux or macOS  
- Bash  
- CTFFIND3  
- `.mrc` micrographs  

CTFFIND3 download:  
https://grigoriefflab.umassmed.edu/ctf_estimation_ctffind_ctftilt

---

## Usage

### 1. Prepare

Place all scripts in your working directory and make them executable:

```bash
chmod +x run_ctffind3.sh run_ctf_fit.sh
```

---

### 2. Run CTFFIND3

Process all micrographs matching:

```
*patch_aligned.mrc
```

```bash
./run_ctffind3.sh
```

Outputs:
- `*_ctf.mrc` (power spectrum)
- `*_ctffind3.log` (CTF fitting log)

Existing log files are skipped automatically.

---

### 3. Calculate CTF fit resolution

```bash
./run_ctf_fit.sh
```

Parses CTFFIND3 power spectrum and log files to compute fitted CTF resolution.

---

## Notes

- Edit CPU settings and CTFFIND3 path in `run_ctffind3.sh`
- Designed for batch processing on workstations or HPC systems
