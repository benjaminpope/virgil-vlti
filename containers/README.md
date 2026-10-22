# ESO pipeline containers

`eso-pipelines.def` is ESO's Apptainer definition file for the VLT pipelines
(<https://eso.org/sci/software/apptainer/eso-pipelines.def>, fetched 2026-10-04:
Fedora plus ESO's RPM repository), with two changes marked in the file. The Fedora release is
a build argument, defaulting to 43: ESO's file pins 41, which reached end of life on
2025-12-15, and ESO's repository also has the pipelines for 42 and 43. EDPS, ESO's
command-line workflow runner, is added so that reductions can run without the
EsoReflex GUI.

The image is built by the `eso-pipelines-sif` GitHub Actions workflow, because
unprivileged builds on some HPC clusters map a single user ID, which is too few for `rpm`
(both attempts failed on 2026-10-04). The workflow runs on any pull request
that changes this directory (with the GRAVITY pipeline), or by hand:

```bash
gh workflow run eso-pipelines-sif.yml -R benjaminpope/virgil-vlti -f pipelines=gravity
```

Each run uploads the `.sif`, its SHA-256, the installed RPM versions and the list
of esorex recipes as an artifact kept for 30 days:

```bash
gh run download <run-id> -R benjaminpope/virgil-vlti -n eso-pipelines-gravity -D eso-pipelines-sif
```

On an HPC cluster, put the image in shared storage. It runs with
`apptainer exec -B <data dir> <image>.sif esorex ...` after `module load apptainer`.
The first image (run 37171175989, from the PR that added this) is on Fedora 41.
