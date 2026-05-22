# Reproducible demo GIF for the README.
#
# Requires `vhs` (https://github.com/charmbracelet/vhs):
#   winget install charmbracelet.vhs            # Windows
#   brew install vhs                            # macOS
#   go install github.com/charmbracelet/vhs@latest
#
# Then from the repo root:
#   vhs docs/demo.tape
# Output: docs/demo.gif

Output docs/demo.gif

Set Shell "bash"
Set FontSize 16
Set Width 1280
Set Height 720
Set Theme "Dracula"
Set TypingSpeed 40ms
Set PlaybackSpeed 1.0

# Peek at the prescription
Type "cat examples/singlet_bk7.yaml"
Enter
Sleep 2500ms

Type "clear"
Enter
Sleep 200ms

# Headline numbers at 5° off-axis
Type "analyze info examples/singlet_bk7.yaml --field-angle-deg 5 --rings 4"
Enter
Sleep 3500ms

# The differentiator: physics self-checks
Type "analyze validate"
Enter
Sleep 3500ms

# And the 4-panel report drops on disk
Type "analyze report examples/singlet_bk7.yaml -o cooke.pdf"
Enter
Sleep 2500ms
