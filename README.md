# iCPSS Risk Score Calculator

Client-side calculator for the **International CMML Prognostic Scoring System (iCPSS)** risk score in Chronic Myelomonocytic Leukemia (CMML). Runs entirely in the browser — no data leaves your device.

## Features

- **Single Patient** — step-by-step form with instant scoring and outcome estimates
- **Cohort Scoring** — download an Excel template, fill it locally, re-upload for batch computation
- Confidence estimate when gene mutation status or karyotype is unknown
- Optional demographic adjustment (age/sex) for patients aged 65–90

## Usage

**Online:** [icpss-risk.com](https://www.icpss-risk.com/)

**Run locally:** download or clone this repository and open `index.html` in any browser. No server or build step required.

```bash
git clone https://github.com/lucalanino/icpss-webapp.git
open icpss-webapp/index.html
```

## Privacy

All computation runs locally in the browser. No patient data is transmitted, stored, or uploaded.

## Reference

Lanino L, Hunter AM, et al. Molecular-Based Ecosystem to Improve Personalized Medicine in Chronic Myelomonocytic Leukemia. *In press*.

## Acknowledgments

- Built with [Claude](https://claude.ai)
- Batch processing uses [SheetJS Community Edition](https://sheetjs.com) (Apache 2.0)
