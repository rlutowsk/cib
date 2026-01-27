#
# cib: Cofinite Integral Braces in GAP
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "cib",
Subtitle := "Cofinite Integral Braces in GAP",
Version := "0.1",
Date := "27/01/2026", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Rafał",
    LastName := "Lutowski",
    WWWHome := "mat.ug.edu.pl/~rlutowsk",
    Email := "rafal.lutowski@ug.edu.pl",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := "80-309",
    Place := "Gdańsk",
    Institution := "University of Gdańsk",
  ),
],

SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/rlutowsk/cib",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://rlutowsk.github.io/cib/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "cib",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Cofinite Integral Braces in GAP",
),

Dependencies := rec(
  GAP := ">= 4.13",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
  if not IsKernelExtensionAvailable("cib") then
    LogPackageLoadingMessage(PACKAGE_WARNING,
                             "failed to load kernel module of package cib");
    return false;
  fi;
  return true;
end,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


