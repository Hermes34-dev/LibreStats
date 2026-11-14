# Releasing LibreStats

Installed copies of LibreStats check for updates at

    https://hermes34-dev.github.io/LibreStats/update.xml

(set in `librestats/src/description.xml`). Every installed copy keeps asking that exact
address, so **never move or rename it**: don't rename the repository, the GitHub account
or the Pages folder. The extension identifier `org.librestats` must never change either:
LibreOffice treats a different identifier as a different extension.

## One-time setup: GitHub Pages

Repository **Settings → Pages → Build and deployment**: Source *Deploy from a branch*,
branch `main`, folder `/docs`. Pages then serves `docs/update.xml` at the address above.

## Each release

1. **Set the version** in `librestats/src/description.xml` (`<version value="X.Y.Z"/>`).
   It is the only place: `build.sh` reads it from there. Also update the version line at
   the top of `documentation.md` and the status line in `README.md`.
2. **Build:** `cd librestats && ./build.sh`. This writes to `librestats/build/`:
   - `librestats.oxt`: the extension
   - `update.xml`: the update feed for this version
   - `SHA256SUMS`: checksum of the `.oxt`
3. **Test** the built `.oxt` (unit tests, and `tests/test_extension.py` after installing it).
4. **Create the GitHub release** with tag `vX.Y.Z` (the `v` matters: the feed links to
   `releases/download/vX.Y.Z/librestats.oxt`). Attach `librestats.oxt` with exactly that
   file name, and `SHA256SUMS`. The release page is also the release notes that the
   Extension Manager links to, so describe the changes there.
5. **Only after the release is published**, copy `librestats/build/update.xml` to
   `docs/update.xml`, commit and push. If the feed goes out first, anyone who checks for
   updates in between gets a broken download.

Check: `https://hermes34-dev.github.io/LibreStats/update.xml` shows the new version
(Pages can take a minute or two), and its download link works.

## How users get the update

In LibreOffice: **Tools → Extension Manager → Check for Updates**. LibreOffice reads the
feed, offers the new version if it is newer than the installed one and the user's
LibreOffice meets the minimum version in the feed, then downloads and installs it.
