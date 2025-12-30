#!/usr/bin/env python3
"""Monitor MEDIUM test progress - memory and status"""

import time
import re
from pathlib import Path

def monitor():
    log_file = Path(r"E:\Projects\halogenator\test_medium_fixed_output.log")

    print("="*80)
    print("MEDIUM TEST MONITOR - Real-time Progress")
    print("="*80)
    print()
    print("Monitoring: test_medium_fixed_output.log")
    print("Press Ctrl+C to stop monitoring")
    print()

    last_size = 0
    last_mem_values = []

    while True:
        try:
            if not log_file.exists():
                print("Waiting for log file...")
                time.sleep(5)
                continue

            with open(log_file) as f:
                content = f.read()

            current_size = len(content)
            if current_size == last_size:
                print(".", end="", flush=True)
                time.sleep(5)
                continue

            last_size = current_size

            # Extract latest memory values
            mem_values = re.findall(r'Pre-flush memory: ([0-9.]+)%', content)
            if mem_values:
                mem_floats = [float(m) for m in mem_values]
                current_max = max(mem_floats[-100:]) if len(mem_floats) > 100 else max(mem_floats)
                current_avg = sum(mem_floats[-100:]) / len(mem_floats[-100:]) if len(mem_floats) > 100 else sum(mem_floats) / len(mem_floats)

                print(f"\rMemory: Current avg={current_avg:.1f}%, max={current_max:.1f}%, flushes={len(mem_floats)}", end="", flush=True)

                # Check for critical warnings
                critical_count = content.count('Critical system memory')
                if critical_count > len(last_mem_values):
                    print(f"\n[WARNING] Critical memory warning detected! Count: {critical_count}")

                last_mem_values = mem_values

            # Check if complete
            if 'TEST COMPLETE' in content or 'TRANSFORMATION COMPLETE' in content:
                print("\n\n" + "="*80)
                print("TEST FINISHED!")
                print("="*80)

                # Final summary
                if mem_values:
                    mem_floats = [float(m) for m in mem_values]
                    print(f"\nFinal Memory Statistics:")
                    print(f"  Peak: {max(mem_floats):.1f}%")
                    print(f"  Average: {sum(mem_floats)/len(mem_floats):.1f}%")
                    print(f"  Total flushes: {len(mem_floats)}")

                critical_count = content.count('Critical system memory')
                print(f"  Critical warnings: {critical_count}")

                print(f"\nFull results in: {log_file}")
                break

            time.sleep(5)

        except KeyboardInterrupt:
            print("\n\nMonitoring stopped by user")
            break
        except Exception as e:
            print(f"\nError: {e}")
            time.sleep(5)

if __name__ == '__main__':
    monitor()
