import os.path
import pathlib
import pytest
import subprocess
import tempfile
from ccbr_tools.shell import shell_run


def extract_command_line(output):
    """Extract the nextflow command line from sinclair run output."""
    cmd_line = None

    # Nextflow 26+ format: ccbr_tools echoes the command as 'nextflow run ...'
    for line in output.split("\n"):
        stripped = line.strip()
        if stripped.startswith("nextflow run") and cmd_line is None:
            cmd_line = stripped

    # Fallback for older Nextflow format: 'cmd line: nextflow run ...'
    if cmd_line is None:
        for line in output.split("\n"):
            if ":" in line and line.split(":")[0].strip() == "cmd line":
                cmd_line = line.split(":", 1)[1].strip()
                break

    if cmd_line is None:
        raise ValueError(f"Could not extract command line from output:\n{output}")

    return cmd_line


def test_help():
    output = subprocess.run(
        "./bin/sinclair --help", capture_output=True, shell=True, text=True
    ).stdout
    assert "Usage: sinclair [OPTIONS]" in output


def test_version():
    output = subprocess.run(
        "./bin/sinclair --version", capture_output=True, shell=True, text=True
    ).stdout
    assert "sinclair, version" in output


def test_citation():
    output = subprocess.run(
        "./bin/sinclair --citation", capture_output=True, shell=True, text=True
    ).stdout
    assert "title = {SINCLAIR" in output


def test_preview():
    with tempfile.TemporaryDirectory() as tmp_dir:
        output = shell_run(
            f"./bin/sinclair init --output {tmp_dir} && ./bin/sinclair run --output {tmp_dir} --mode local -preview -profile ci_stub"
        )
    cmd_line = extract_command_line(output)
    assert all(["-preview" in cmd_line, "-resume" in cmd_line])


def test_forceall():
    output = subprocess.run(
        "./bin/sinclair run --forceall -preview -profile ci_stub --mode local",
        capture_output=True,
        shell=True,
        text=True,
        check=True,
    ).stdout
    cmd_line = extract_command_line(output)
    assert "-preview" in cmd_line and "-resume" not in cmd_line


def test_init():
    with tempfile.TemporaryDirectory() as tmp_dir:
        shell_run(f"./bin/sinclair init --output {tmp_dir}")
        outdir = pathlib.Path(tmp_dir)
        assertions = [(outdir / "nextflow.config").exists(), (outdir / "log").exists()]
    assert all(assertions)


def test_init_default():
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        shell_run(f"{cwd}/bin/sinclair init")
        outdir = pathlib.Path(tmp_dir)
        assertions = [(outdir / "nextflow.config").exists(), (outdir / "log").exists()]

    os.chdir(cwd)
    assert all(assertions)


def test_run_no_init():
    with pytest.raises(Exception):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = shell_run(
                f"./bin/sinclair run --output {tmp_dir} --mode local",
                check=True,
                capture_output=True,
            )
            assertions = ["Hint: you must initialize the output directory" in output]
            assert all(assertions)
