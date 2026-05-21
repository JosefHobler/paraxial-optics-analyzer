"""Pin matplotlib to a non-interactive backend before pyplot imports anywhere.

Without this, Windows runs default to TkAgg, which intermittently fails to
locate init.tcl in headless / pytest contexts. The report writes PDFs/PNGs to
disk — no GUI needed.
"""
import matplotlib

matplotlib.use("Agg")
