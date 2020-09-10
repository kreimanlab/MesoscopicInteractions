#!/bin/bash

# Download
echo rsync -avP ../downloads.tar.gz ./
rsync -avP ../downloads.tar.gz ./
echo rm -rf downloads
rm -rf downloads
echo tar -xzf downloads.tar.gz
echo "  (*) Decompressing.."
tar -xzf downloads.tar.gz
#echo rm downloads.tar.gz
#rm downloads.tar.gz

# Extract example .txt
echo cp ./data/txt/example/archive/_Export-example_1.txt.xz ./data/txt/example/
cp ./data/txt/example/archive/_Export-example_1.txt.xz ./data/txt/example/
echo rm ./data/txt/example/_Export-example_1.txt
rm ./data/txt/example/_Export-example_1.txt
echo unxz ./data/txt/example/_Export-example_1.txt.xz
unxz ./data/txt/example/_Export-example_1.txt.xz

echo "[!] Please check contents of ./downloads/ folder:"
echo "  - Install freesurfer if it is not already on your system."
echo "      - Configure \$SUBJECTS_DIR to ./coregistration/"
echo "  - Install iElVis if it is not already on your system."
