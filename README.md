# ULM Threshold Analysis

This repository contains code and scripts used for analyzing the effect of detection thresholds—specifically false positives and false negatives—on **Super-Resolution Ultrasound (SR-US)** and **Ultrasound Localization Microscopy (ULM)** image quality.

The focus of this work is to quantitatively study how microbubble detection errors propagate into super-resolution density maps and affect commonly used image quality metrics.

---

## 📦 Data Availability

The datasets used are from UltraSR challenge and can be downloaded from **Zenodo**:

🔗 **https://zenodo.org/records/7271766**

Please download the data separately and place it in the appropriate directory structure as expected by the scripts in this repository.

---

## 📊 What This Repository Includes

- Scripts for modifying microbubble detection outputs by:
  - Adding controlled numbers of **false positives**
  - Introducing **false negatives** at specified rates
- Generation of **super-resolution density maps**
- Quantitative evaluation using image quality metrics such as:
  - Structural Similarity Index (SSIM)
  - Dice coefficient
  - Peak Signal-to-Noise Ratio (PSNR)
- Analysis of dense vs. sparse regions in SR maps

The code is intended for **research and reproducibility purposes** and may require adaptation for different datasets or experimental setups.

---
## 📚 Citation

If you use this code or data in your work, **please cite the following paper for now**:

> **ArXiv preprint**  
> https://arxiv.org/abs/2411.07426

⚠️ **Note:**  
This work has been **accepted to SPIE Medical Imaging 2026**.  
The citation information will be updated once the final published version becomes available.

---

## 📬 Contact

For questions, issues, or collaboration inquiries, please open an issue on GitHub or contact the authors directly.

---

## 📄 License

This repository is intended for academic and research use.  
Please check the Zenodo record for data-specific licensing information.

---

