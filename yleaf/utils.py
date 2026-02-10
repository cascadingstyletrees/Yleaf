import logging
import subprocess
import shlex
from typing import IO, Any

from yleaf.exceptions import ExternalCommandError

LOG = logging.getLogger("yleaf_logger")

def run_command(
    command: str | list[str],
    stdout: int | IO[Any] | None = None,
    shell: bool = False,
    check: bool = True,
    capture_output: bool = False,
    text: bool = True
) -> subprocess.CompletedProcess:
    """
    Runs a shell command using subprocess.run.

    Args:
        command: Command to run. Can be a string or a list of arguments.
        stdout: File object to write stdout to. If None, stdout is captured if capture_output is True.
        shell: Whether to run the command in a shell.
        check: Whether to raise an exception if the command fails.
        capture_output: Whether to capture stdout and stderr (if not redirecting to file).
        text: Whether to treat output as text (decode bytes).

    Returns:
        The CompletedProcess object.

    Raises:
        ExternalCommandError: If the command fails and check is True.
    """

    # Ensure command is in correct format for shell argument
    if shell and isinstance(command, list):
        # Convert list to string for shell execution
        # Use shlex.join if available (Python 3.8+), otherwise manual join
        # But simple join " ".join(command) is what Yleaf did before.
        # shlex.join handles quoting which is safer.
        cmd_for_run = " ".join(command)
        cmd_log_str = cmd_for_run
    else:
        cmd_for_run = command
        if isinstance(command, list):
            cmd_log_str = " ".join(command)
        else:
            cmd_log_str = command

    LOG.debug(f"Running command: {cmd_log_str}")

    try:
        # Determine stderr handling: always capture for error reporting
        stderr_dest = subprocess.PIPE

        # Determine stdout handling
        stdout_dest = stdout
        if stdout is None and capture_output:
            stdout_dest = subprocess.PIPE

        result = subprocess.run(
            cmd_for_run,
            shell=shell,
            check=False, # Check manually
            stdout=stdout_dest,
            stderr=stderr_dest,
            text=text
        )

        if check and result.returncode != 0:
            stderr_output = result.stderr if result.stderr else "No stderr captured."
            LOG.error(f"Command failed with return code {result.returncode}")
            # Log stderr if it wasn't empty
            if stderr_output.strip():
                LOG.error(f"Stderr: {stderr_output.strip()}")

            raise ExternalCommandError(cmd_log_str, result.returncode, str(stderr_output))

        return result

    except FileNotFoundError as e:
        # Command not found (only happens if shell=False)
        LOG.error(f"Command not found: {e}")
        raise ExternalCommandError(cmd_log_str, 127, f"Command not found: {e}")
    except Exception as e:
        # Catch unexpected errors
        if isinstance(e, ExternalCommandError):
            raise
        LOG.error(f"Unexpected error running command: {e}")
        raise ExternalCommandError(cmd_log_str, -1, str(e))
