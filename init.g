#
# cib: Cofinite Integral Braces in GAP
#
# Reading the declaration part of the package.
#
if not LoadKernelExtension("cib") then
  Error("failed to load kernel module of package cib");
fi;
ReadPackage( "cib", "gap/cib.gd");
