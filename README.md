# fuzzig

Fuzzy finder algorithms in Zig based on the [Smith-Waterman algorithm](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm), inspired by [fzf](https://github.com/junegunn/fzf).

## Example

```zig
const std = @import("std");
const fuzzig = @import("fuzzig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // The default implementation needs to know the maximum haystack and needle
    // sizes so that it can allocate all memory contiguously. Smaller haystack
    // or needle strings will use subsets of the allocated memory.
    var searcher = try fuzzig.Ascii.init(
        allocator,
        128, // haystack max size
        32, // needle max size
        .{ .case_sensitive = false },
    );
    defer searcher.deinit();

    const haystack = "Hello World";
    const needle = "world";

    const score = searcher.score(haystack, needle);
    std.debug.print("Score: {d}\n", .{score.?});

    // Get a traceback of the character positions that were matched
    const match = searcher.scoreMatches(haystack, needle);
    std.debug.print(
        "Score with traceback: {d} {any}\n",
        .{ match.score.?, match.matches },
    );
}
```

Output

```
Score: 104
Score with traceback: 104 { 6, 7, 8, 9, 10 }
```

See the [`AsciiOptions` struct](https://github.com/fjebaker/fuzzig/blob/a78afddec30b547643604aafaee202db6fc878f1/src/root.zig#L460-L466) for a list of available options.

## Design

The module defines an `Algorithm` generic type, which accepts the element type of the array to be fuzzy searched, the score type and values, and an algorithm implementation. The implementation must define an `eqlFunc`, a `scoreFunc` and a `bonusFunc` used to test for equality between tokens, for determining the score of two matching tokens, and for determining any in-places bonuses respectively.

This package currently only implements an [ASCII fuzzy finder](https://github.com/fjebaker/fuzzig/blob/a78afddec30b547643604aafaee202db6fc878f1/src/root.zig#L457-L511), but can be used as reference to extend to e.g. UTF8 or other character encodings.

## Usage

The library was written with Zig 0.12.0-dev.3541+05b185811, but there is likely a lot of flexibility with versioning.

To use in a Zig project, add it to your `build.zig.zon`

```zig
    // ...
    .dependencies = .{
        .fuzzig = .{
            .url = "https://github.com/fjebaker/fuzzig/archive/main.tar.gz",
            .hash = "" // get with `zig fetch`
        },
    },
    // ...
```

Then add the module to your build step in `build.zig`:

```zig
    // ...
    const fuzzig = b.dependency("fuzzig", .{}).module("fuzzig");

    my_exe_or_lib.root_module.addImport("fuzzig", fuzzig);
    // ...
```
