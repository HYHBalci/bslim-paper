$ErrorActionPreference = 'Stop'
Push-Location -LiteralPath $PSScriptRoot
try {
    $documents = @('methodology-8.3.5', 'methodology-8.3.5-supplement')
    foreach ($document in $documents) {
        & pdflatex -interaction=nonstopmode -halt-on-error "$document.tex"
        if ($LASTEXITCODE -ne 0) { throw "LaTeX failed: $document" }
        & bibtex $document
        if ($LASTEXITCODE -ne 0) { throw "BibTeX failed: $document" }
    }
    # Both documents refer to labels in the other document.
    foreach ($pass in 1..3) {
        foreach ($document in $documents) {
            & pdflatex -interaction=nonstopmode -halt-on-error "$document.tex"
            if ($LASTEXITCODE -ne 0) { throw "LaTeX failed: $document" }
        }
    }
} finally {
    Pop-Location
}
