# Code signing policy

The Windows installer of SangerScope is built from the source code in this
repository by the "Windows installer" GitHub Actions workflow
(`.github/workflows/windows-installer.yml`). Only installers produced by that
workflow from the `master` branch are published as releases.

Free code signing provided by [SignPath.io](https://signpath.io), certificate by
[SignPath Foundation](https://signpath.org).

## Team roles

- Committers and reviewers: [alrafaykhan](https://github.com/alrafaykhan)
- Approvers of signing requests: [alrafaykhan](https://github.com/alrafaykhan)

## Privacy

SangerScope runs on the user's own computer. Trace files and results are not
sent anywhere. The program contacts other servers only in these cases:

- NCBI, when the user searches for or downloads a reference sequence: only the
  search words or the accession are sent.
- GitHub, once on start, to ask whether a newer version exists. This can be
  switched off by setting `SANGERSCOPE_NO_UPDATE_CHECK=1`.
- NCBI Variation Services, when a human variant is shown, to find its position on the
  genome for the Franklin link: only the variant name (transcript and c. change) is sent.
- ClinVar and Franklin, only when the user clicks the link of a variant.
