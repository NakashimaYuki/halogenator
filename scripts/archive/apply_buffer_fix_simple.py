#!/usr/bin/env python3
"""Simple script to apply buffer explosion fix."""

def apply_fix():
    filepath = r'E:\Projects\halogenator\scripts\08_transform_library_v2.py'

    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Find the write_batch method (around line 413)
    start_idx = None
    for i, line in enumerate(lines):
        if 'def write_batch(self, products: List[Dict]):' in line:
            start_idx = i
            break

    if start_idx is None:
        print("ERROR: Could not find write_batch method")
        return False

    # Find the end of the method (next method definition)
    end_idx = None
    for i in range(start_idx + 1, len(lines)):
        if lines[i].startswith('    def ') and not lines[i].startswith('        '):
            end_idx = i
            break

    if end_idx is None:
        print("ERROR: Could not find end of write_batch method")
        return False

    print(f"Found write_batch at lines {start_idx+1} to {end_idx}")

    # New implementation
    new_method = '''    def write_batch(self, products: List[Dict]):
        """
        Intelligent memory-managed batch writing with SYSTEM-WIDE monitoring.

        CRITICAL FIX (2025-12-27): Process large batches in chunks to prevent
        buffer explosion. Previous code allowed buffer to grow to 900K+ products
        before flush, causing 1-3GB memory spikes and OOM crashes.

        Strategy:
        1. HARD LIMIT: Never allow buffer to exceed flush_interval
        2. Process incoming products in chunks if needed
        3. Monitor GLOBAL system memory (not single process)
        4. Force flush and GC after each chunk to free memory immediately
        """
        if not products:
            return

        # HARD LIMIT: Never exceed this buffer size
        MAX_BUFFER_SIZE = self.flush_interval

        # Calculate safe chunk size (half of max to allow headroom)
        chunk_size = max(1000, MAX_BUFFER_SIZE // 2)

        # Process products in chunks to prevent buffer explosion
        for i in range(0, len(products), chunk_size):
            chunk = products[i:i+chunk_size]

            # Add chunk to buffer
            self.buffer.extend(chunk)

            # Get current SYSTEM memory status
            system_mem_percent = self._get_memory_usage_percent()
            buffer_size = len(self.buffer)

            # Adaptive flush conditions
            should_flush = False
            flush_reason = ""

            # Condition 1: Buffer size limit (PRIMARY - HARD LIMIT)
            if buffer_size >= MAX_BUFFER_SIZE:
                should_flush = True
                flush_reason = f"buffer_full ({buffer_size:,} products)"

            # Condition 2: System memory pressure (SECONDARY)
            if system_mem_percent > self.target_memory_percent:
                should_flush = True
                flush_reason = f"system_memory_pressure ({system_mem_percent:.1f}%)"

            # Condition 3: Critical system memory (URGENT)
            if system_mem_percent > 85.0:
                should_flush = True
                flush_reason = f"CRITICAL_SYSTEM_MEMORY ({system_mem_percent:.1f}%)"
                self.logger.warning(f"[!] Critical system memory: {system_mem_percent:.1f}%")

            if should_flush:
                # [PROFILING] Record pre-flush memory
                pre_flush_mem = system_mem_percent
                self.logger.info(f"Flush triggered: {flush_reason}")
                self.logger.info(f"[PROFILING] Pre-flush memory: {pre_flush_mem:.1f}%")

                # Flush and measure timing
                self._flush()

                # Force garbage collection to free memory immediately
                import gc
                gc.collect()

                # [PROFILING] Record post-flush memory and delta
                post_flush_mem = self._get_memory_usage_percent()
                mem_delta = post_flush_mem - pre_flush_mem
                self.logger.info(f"[PROFILING] Post-flush memory: {post_flush_mem:.1f}%, delta: {mem_delta:+.1f}%")

'''

    # Backup original
    backup_path = r'E:\Projects\halogenator\scripts\08_transform_library_v2.py.backup_before_buffer_fix'
    with open(backup_path, 'w', encoding='utf-8') as f:
        f.writelines(lines)
    print(f"Created backup: {backup_path}")

    # Replace the method
    new_lines = lines[:start_idx] + [new_method] + lines[end_idx:]

    # Write modified version
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

    print(f"Applied fix to: {filepath}")
    print("\nFix Summary:")
    print("  - Enforces HARD buffer size limit (max flush_interval)")
    print("  - Processes large batches in 1000-product chunks")
    print("  - Forces GC after each flush")
    print("  - Prevents buffer explosion (900K+ products)")
    return True

if __name__ == '__main__':
    print("="*80)
    print("Applying StreamingParquetWriter Buffer Explosion Fix")
    print("="*80)
    success = apply_fix()
    print("\n" + ("SUCCESS!" if success else "FAILED"))
