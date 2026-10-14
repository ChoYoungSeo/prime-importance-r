"""Interrupt only this script's fresh R subprocess; enforce an 8-second bound."""
import argparse
import json
import os
from pathlib import Path
import selectors
import signal
import subprocess
import time

folder = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--library", type=Path, required=True)
parser.add_argument("--rscript", default="Rscript")
args = parser.parse_args()
env = os.environ.copy()
for name in ("PRIME_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
             "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    env[name] = "1"
started = time.monotonic()
process = subprocess.Popen(
    [args.rscript, "--vanilla", str(folder / "interrupt_worker.R"),
     str(args.library.resolve())],
    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env,
)
lines = []
try:
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    while time.monotonic() - started < 3:
        if selector.select(timeout=0.1):
            line = process.stdout.readline()
            lines.append(line)
            if "READY_FOR_OWN_PROCESS_SIGINT" in line:
                break
            if not line and process.poll() is not None:
                raise RuntimeError("R exited before entering native workload")
    else:
        raise TimeoutError("R did not reach native workload within 3 seconds")
    # Allow the C++ entry point to launch its two candidate workers.
    time.sleep(0.15)
    sent = time.monotonic()
    process.send_signal(signal.SIGINT)
    output, _ = process.communicate(timeout=max(0.1, 8-(sent-started)))
    lines.append(output)
    latency = time.monotonic() - sent
    transcript = "".join(lines)
    result = {
        "subprocess_pid": process.pid,
        "exit_code": process.returncode,
        "signal": "SIGINT sent only to spawned R process",
        "signal_to_exit_seconds": latency,
        "total_seconds": time.monotonic()-started,
        "clean_interrupt_caught": "CAUGHT_INTERRUPT_CLEANLY" in transcript,
    }
    assert process.returncode == 0 and result["clean_interrupt_caught"], transcript
    (folder / "interrupt_result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(transcript, end="")
    print(json.dumps(result, indent=2))
finally:
    if process.poll() is None:
        process.kill()
        process.communicate(timeout=1)
