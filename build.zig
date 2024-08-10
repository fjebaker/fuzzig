const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zg = b.dependency("zg", .{
        .target = target,
        .optimize = optimize,
    });

    _ = b.addModule("fuzzig", .{ .root_source_file = b.path("src/root.zig") });

    const lib = b.addStaticLibrary(.{
        .name = "fuzzig",
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    lib.root_module.addImport("code_point", zg.module("code_point"));
    lib.root_module.addImport("GenCatData", zg.module("GenCatData"));
    lib.root_module.addImport("CaseData", zg.module("CaseData"));
    lib.root_module.addImport("Normalize", zg.module("Normalize"));
    lib.root_module.addImport("CaseFold", zg.module("CaseFold"));

    b.installArtifact(lib);

    const lib_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    lib_unit_tests.root_module.addImport("code_point", zg.module("code_point"));
    lib_unit_tests.root_module.addImport("GenCatData", zg.module("GenCatData"));
    lib_unit_tests.root_module.addImport("CaseData", zg.module("CaseData"));
    lib_unit_tests.root_module.addImport("Normalize", zg.module("Normalize"));
    lib_unit_tests.root_module.addImport("CaseFold", zg.module("CaseFold"));

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    const exe = b.addExecutable(.{
        .name = "fuzzig",
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    const run_cmd = b.addRunArtifact(exe);
    const benchmark_step = b.step("benchmark", "Run benchmarks.");
    benchmark_step.dependOn(&run_cmd.step);
}
