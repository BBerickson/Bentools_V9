# Small shared UI helpers.

# DT column renderer: show the full cell text as a hover tooltip and truncate the
# visible text. `threshold` = only add the tooltip/truncation when the value is
# longer than this many characters; `len` = number of characters shown before "…".
# Replaces a JS(...) block that was duplicated across ~12 datatable() calls.
trunc_render <- function(threshold, len) {
  htmlwidgets::JS(
    "function(data, type, row, meta) {",
    sprintf("return type === 'display' && data.length > %d ?", threshold),
    sprintf("'<span title=\"' + data + '\">' + data.substr(0, %d) + '...</span>' : data;", len),
    "}"
  )
}
