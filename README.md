# iCPSS Risk Score Calculator

<<<<<<< HEAD
Client-side calculator for the **Integrated Clinical and Prognostic Scoring System (iCPSS)** risk score in Chronic Myelomonocytic Leukemia (CMML). Runs entirely in the browser — no data is uploaded or transmitted.
=======
Client-side calculator for the **International Clinical and Prognostic Scoring System (iCPSS)** risk score in Chronic Myelomonocytic Leukemia (CMML). Runs entirely in the browser — no data is uploaded or transmitted.
>>>>>>> 6628b67 (Update README)

## Features

- **Single Patient** — interactive form with instant scoring and outcome estimates
- **Cohort Scoring** — download an Excel template, fill it locally, re-upload for batch computation
- Confidence estimate when gene mutation status or karyotype is unknown
- Optional demographic adjustment (age/sex) for patients aged 65–90

## Usage

Open `index.html` in any modern browser. No server or build step required.

## Privacy

No external network requests. No analytics. No data upload of any kind.

## Reference

*[Citation in press]*

## Acknowledgments

Batch processing uses [SheetJS Community Edition](https://sheetjs.com) (Apache 2.0).