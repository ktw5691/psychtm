# psychtm

## Troubleshooting `C++` on Mac OS

**I am not responsible for what happens if you try the following:**

Using `R` packages that rely on `C++` code can be tricky on Mac OS. Common
errors include warnings during package installation that certain files (e.g.,
`math.h`) are unavailable. This can generally be resolved by executing the
following code in the terminal:

```{bash}
xcode-select --install
installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
```

The first line `xcode-select --install` will install the Command Line Tools
for XCode if you do not already have them installed which will be necessary for
many uses of `C++`.

The second line `installer -pkg ...` installs a variety of `C++` header files
(including `math.h`) that are needed for compiling many `C++` programs and are
not always available by default on Mac OS. You may need to execute this line as:

```{bash}
sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
```

to give sufficient permissions for this command.
