# Output files
x_file="x.dat"
y_file="y.dat"
z_file="z.dat"
output="2d-halo-ener.dat"

# Clear temp files if they exist
> "$x_file"
> "$y_file"
> "$z_file"

for np in ARG93 ARG97 ASN171 ASN362 ASN89 FRA440 GLU363 GLY173 GLY95 HIS348 HIS359 HIS96 LEU188 LYS356 LYS357 LYS92 MET94 THR168 TRP98 TYR360 TYR361 VAL347 VAL358
do
  # Get column 3 or blank
  val_x=$(awk '{print ($3 != "" ? $3*10 : "")}' "${np}.dat" 2>/dev/null | head -1)
  echo "${val_x}" >> "$x_file"

  # Get column 4 or blank
  val_y=$(awk '{print ($4 != "" ? $4 : "")}' "${np}.dat" 2>/dev/null | head -1)
  echo "${val_y}" >> "$y_file"

  # Get max of column 15 or blank
  if [ -f "${np}-2PT-E.dat" ]; then
    val_z=$(awk '{if ($15 != "") print $15}' "${np}-2PT-E.dat" | sort -n | tail -1)
    echo "${val_z}" >> "$z_file"
  else
    echo "" >> "$z_file"
  fi
done

# Combine and clean
paste "$x_file" "$y_file" "$z_file" > "$output"
rm "$x_file" "$y_file" "$z_file"

