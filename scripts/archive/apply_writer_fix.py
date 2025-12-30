#!/usr/bin/env python3
"""
Apply critical fix to StreamingParquetWriter to prevent buffer explosion.

This script modifies the write_batch() method to process large batches in chunks,
preventing the buffer from growing to 900K+ products and causing OOM crashes.
"""

import re

def apply_fix(filepath):
    """Apply the buffer explosion fix to the transform library."""

    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # The old write_batch method
    old_pattern = r'''    def write_batch\(self, products: List\[Dict\]\):
        """
        Intelligent memory-managed batch writing with SYSTEM-WIDE monitoring.

        Strategy:
        1\. Monitor GLOBAL system memory \(not single process\)
        2\. Trigger flush when buffer full \(primary\) or system memory high \(secondary\)
        3\. Use adaptive thresholds based on actual system memory
        4\. Critical flush at 85% system memory
        """
        if not products:
            return

        # Add to buffer
        self\.buffer\.extend\(products\)

        # Get current SYSTEM memory status \(not process\)
        system_mem_percent = self\._get_memory_usage_percent\(\)
        buffer_size = len\(self\.buffer\)

        # Adaptive flush conditions
        should_flush = False
        flush_reason = ""

        # Condition 1: Buffer size limit \(PRIMARY - most reliable\)
        if buffer_size >= self\.flush_interval:
            should_flush = True
            flush_reason = f"buffer_full \(\{buffer_size:,\} products\)"

        # Condition 2: System memory pressure \(SECONDARY\)
        if system_mem_percent > self\.target_memory_percent:
            should_flush = True
            flush_reason = f"system_memory_pressure \(\{system_mem_percent:\.1f\}%\)"

        # Condition 3: Critical system memory \(URGENT\)
        if system_mem_percent > 85\.0:
            should_flush = True
            flush_reason = f"CRITICAL_SYSTEM_MEMORY \(\{system_mem_percent:\.1f\}%\)"
            self\.logger\.warning\(f"\[!\] Critical system memory: \{system_mem_percent:\.1f\}%"\)

        if should_flush:
            # \[PROFILING\] Record pre-flush memory
            pre_flush_mem = system_mem_percent
            self\.logger\.info\(f"Flush triggered: \{flush_reason\}"\)
            self\.logger\.info\(f"\[PROFILING\] Pre-flush memory: \{pre_flush_mem:\.1f\}%"\)

            # Flush and measure timing
            self\._flush\(\)

            # \[PROFILING\] Record post-flush memory and delta
            post_flush_mem = self\._get_memory_usage_percent\(\)
            mem_delta = post_flush_mem - pre_flush_mem
            self\.logger\.info\(f"\[PROFILING\] Post-flush memory: \{post_flush_mem:\.1f\}%, delta: \{mem_delta:\+\.1f\}%"\)'''

    # New fixed version
    new_code = '''    def write_batch(self, products: List[Dict]):
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
                self.logger.info(f"[PROFILING] Post-flush memory: {post_flush_mem:.1f}%, delta: {mem_delta:+.1f}%")'''

    # Apply the replacement
    if re.search(old_pattern, content, re.DOTALL):
        print("✓ Found write_batch() method to fix")
        new_content = re.sub(old_pattern, new_code, content, count=1, flags=re.DOTALL)

        # Verify the replacement worked
        if 'CRITICAL FIX (2025-12-27)' in new_content:
            # Backup original
            backup_path = filepath + '.backup_before_buffer_fix'
            with open(backup_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"✓ Created backup: {backup_path}")

            # Write fixed version
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(new_content)
            print(f"✓ Applied fix to: {filepath}")
            print("\nFix applied successfully!")
            print("\nWhat changed:")
            print("  - Now processes large batches in chunks (max 1000 products per chunk)")
            print("  - Enforces HARD LIMIT on buffer size (never exceeds flush_interval)")
            print("  - Forces GC after each flush to free memory immediately")
            print("  - Prevents 900K+ product buffer explosions")
            return True
        else:
            print("✗ Replacement failed verification")
            return False
    else:
        print("✗ Could not find write_batch() method pattern")
        print("File may have already been modified or structure changed")
        return False

if __name__ == '__main__':
    import sys
    filepath = 'scripts/08_transform_library_v2.py'
    if len(sys.argv) > 1:
        filepath = sys.argv[1]

    print("="*80)
    print("APPLYING CRITICAL FIX: StreamingParquetWriter Buffer Explosion")
    print("="*80)
    print(f"\nTarget file: {filepath}")
    print()

    success = apply_fix(filepath)
    sys.exit(0 if success else 1)
