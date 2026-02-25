## 2024-10-25 - CLI Download Progress
**Learning:** Users lack feedback during long-running blocking operations like downloading large reference files, leading to uncertainty if the process has hung.
**Action:** Always include progress indicators for potentially long network or I/O operations, especially in CLI tools where silence is ambiguous.
