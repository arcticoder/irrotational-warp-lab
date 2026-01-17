#!/usr/bin/env python3
"""
GPU availability check for irrotational-warp CLI tools.

Verifies CuPy + CUDA setup for GPU-accelerated reproduction runs.
"""

import subprocess
import sys


def check_gpu():
    """Check GPU availability and provide setup guidance."""
    print("=" * 60)
    print("GPU ACCELERATION CHECK")
    print("=" * 60)

    # Check nvidia-smi
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            gpu_info = result.stdout.strip()
            print(f"✅ GPU detected: {gpu_info}")
        else:
            print("❌ nvidia-smi failed - GPU may not be accessible")
            return False
    except FileNotFoundError:
        print("❌ nvidia-smi not found - NVIDIA drivers not installed")
        print("\nWSL2 GPU Setup:")
        print("1. Install NVIDIA GPU drivers on Windows (version 455.41+)")
        print("2. In WSL2, verify GPU: nvidia-smi")
        return False
    except Exception as e:
        print(f"❌ nvidia-smi error: {e}")
        return False

    # Check CuPy
    try:
        import cupy as cp

        print(f"✅ CuPy installed: {cp.__version__}")

        # Test basic GPU operation
        device_count = cp.cuda.runtime.getDeviceCount()
        print(f"✅ CUDA devices available: {device_count}")

        # Quick computation test
        a = cp.array([1, 2, 3, 4, 5])
        b = cp.sum(a)
        result = b.get()
        print(f"✅ GPU computation test passed (sum={result})")

        # Memory info
        mempool = cp.get_default_memory_pool()
        print(f"✅ GPU memory pool initialized")

        print("\n" + "=" * 60)
        print("GPU ACCELERATION READY")
        print("=" * 60)
        print("Use --backend cupy in CLI commands for GPU acceleration")
        return True

    except ImportError:
        print("❌ CuPy not installed")
        print("\nInstall CuPy for CUDA 12.x:")
        print("  pip install cupy-cuda12x")
        print("\nFor other CUDA versions, see:")
        print("  https://docs.cupy.dev/en/stable/install.html")
        return False
    except Exception as e:
        print(f"❌ CuPy error: {e}")
        print("\nPossible issues:")
        print("- CUDA toolkit version mismatch")
        print("- WSL2 GPU passthrough not configured")
        print("- Insufficient GPU memory")
        return False


if __name__ == "__main__":
    success = check_gpu()
    sys.exit(0 if success else 1)
