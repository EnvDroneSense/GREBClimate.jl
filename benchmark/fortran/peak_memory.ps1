# Run one process and print "<peak memory bytes> <exit code>".
# Used by run_fortran_comparison.jl's `memory` mode.
param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$WorkDir,
    [Parameter(Mandatory = $true)][string]$OutFile,
    [string]$PrependPath = "",
    [string]$OmpThreads = "",
    [string]$ArgsFile = ""
)

$ErrorActionPreference = "Stop"

Add-Type -TypeDefinition @"
using System;
using System.Runtime.InteropServices;
public struct PROCESS_MEMORY_COUNTERS {
    public uint cb;
    public uint PageFaultCount;
    public UIntPtr PeakWorkingSetSize;
    public UIntPtr WorkingSetSize;
    public UIntPtr QuotaPeakPagedPoolUsage;
    public UIntPtr QuotaPagedPoolUsage;
    public UIntPtr QuotaPeakNonPagedPoolUsage;
    public UIntPtr QuotaNonPagedPoolUsage;
    public UIntPtr PagefileUsage;
    public UIntPtr PeakPagefileUsage;
}
public static class Psapi {
    [DllImport("psapi.dll", SetLastError=true)]
    public static extern bool GetProcessMemoryInfo(IntPtr hProcess,
        out PROCESS_MEMORY_COUNTERS counters, uint size);
}
"@ -ErrorAction SilentlyContinue

$ExeArgs = @()
if ($ArgsFile -ne "" -and (Test-Path $ArgsFile)) {
    $ExeArgs = @(Get-Content -LiteralPath $ArgsFile | Where-Object { $_ -ne "" })
    # Quote arguments containing spaces.
    $ExeArgs = @($ExeArgs | ForEach-Object {
        if ($_ -match '\s' -and -not $_.StartsWith('"')) { '"' + $_ + '"' } else { $_ }
    })
}

# PrependPath holds the mingw runtime DLLs the Fortran needs.
if ($PrependPath -ne "") { $env:PATH = "$PrependPath;$env:PATH" }
if ($OmpThreads -ne "") { $env:OMP_NUM_THREADS = $OmpThreads }

$opts = @{
    FilePath               = $Exe
    WorkingDirectory       = $WorkDir
    PassThru               = $true
    NoNewWindow            = $true
    RedirectStandardOutput = $OutFile
    RedirectStandardError  = "$OutFile.err"
}
if ($ExeArgs.Count -gt 0) { $opts['ArgumentList'] = $ExeArgs }

$p = Start-Process @opts

# Poll as a fallback.
$polled = 0
while (-not $p.HasExited) {
    try {
        $p.Refresh()
        if ($p.PeakWorkingSet64 -gt $polled) { $polled = $p.PeakWorkingSet64 }
    } catch { }
    Start-Sleep -Milliseconds 250
}
$p.WaitForExit()

# Lifetime peak from the kernel, read after exit.
$peak = 0
try {
    $counters = New-Object PROCESS_MEMORY_COUNTERS
    $size = [System.Runtime.InteropServices.Marshal]::SizeOf($counters)
    if ([Psapi]::GetProcessMemoryInfo($p.Handle, [ref]$counters, $size)) {
        $peak = [uint64]$counters.PeakWorkingSetSize
    }
} catch { }

if ($polled -gt $peak) { $peak = $polled }

[Console]::Out.WriteLine("$peak $($p.ExitCode)")
