using PackageCompiler

create_app(
    ".", 
    "build",
    executables = ["milk" => "julia_main"],
    force = true,
)
