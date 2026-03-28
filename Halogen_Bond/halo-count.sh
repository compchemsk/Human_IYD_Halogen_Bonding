for np in 3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.8 8.0 8.2 8.4 8.6 8.8 9.8 10.8 11.8 12.8 13.8 14.8 15.8 16.8 17.8 18.8
do
    next_np=$(awk -v n=$np 'BEGIN{found=0}
    {if($0 ~ "For trajectory " n) found=1}
    found && match($0, /For trajectory ([0-9]+\.[0-9]+)/, a) {
        if(a[1] > n) {print a[1]; exit}
    }' halo.out)

    if [ -n "$next_np" ]; then
        lines=$(awk "/For trajectory ${np}/{f=1;next} /For trajectory ${next_np}/{f=0} f" halo.out | wc -l)
        echo "Lines between For trajectory $np and For trajectory $next_np: $lines" >> halo-count.dat
    fi
done
cat halo-count.dat | awk '{print($10/500)}' | cat -n > halo-fraction.dat
sed -i 's/0.002/0.0000/g' halo-fraction.dat
echo "    36  0.144" >> halo-fraction.dat
rm halo-count.dat
