"""Refuse to publish unless the tests workflow passed on the commit being released.

Reads GITHUB_TOKEN, GITHUB_REPOSITORY and GITHUB_SHA from the environment. Passes only
if the latest completed run of tests.yml for that commit succeeded, and every job in it
(the Python matrix and the virgil-main job) concluded success.
"""

import json
import os
import sys
import urllib.request

API = "https://api.github.com"
WORKFLOW = "tests.yml"


def get(path):
    request = urllib.request.Request(
        f"{API}{path}",
        headers={
            "Authorization": f"Bearer {os.environ['GITHUB_TOKEN']}",
            "Accept": "application/vnd.github+json",
        },
    )
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def main():
    repo = os.environ["GITHUB_REPOSITORY"]
    sha = os.environ["GITHUB_SHA"]
    url = f"/repos/{repo}/actions/workflows/{WORKFLOW}/runs?head_sha={sha}&per_page=20"
    # Push runs only: a pull_request run reports the PR head but tests the merge with main.
    runs = [
        r for r in get(url)["workflow_runs"] if r["status"] == "completed" and r["event"] == "push"
    ]
    if not runs:
        sys.exit(f"No completed push {WORKFLOW} run for commit {sha}; refusing to publish.")
    run = max(runs, key=lambda r: r["run_started_at"])
    if run["conclusion"] != "success":
        sys.exit(f"{run['html_url']} concluded {run['conclusion']}; refusing to publish.")
    jobs = get(f"/repos/{repo}/actions/runs/{run['id']}/jobs?per_page=100")["jobs"]
    bad = [j["name"] for j in jobs if j["conclusion"] != "success"]
    if not jobs or bad:
        sys.exit(f"Jobs not successful in {run['html_url']}: {bad or 'none'}; not publishing.")
    print(f"{WORKFLOW} passed on {sha} ({len(jobs)} jobs): {run['html_url']}")


if __name__ == "__main__":
    main()
