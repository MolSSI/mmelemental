if [ -f _build ]; then rm -r _build; fi
if [ -f build ]; then rm -r build; fi

make html
cd _build/html

for d in _* ; do
    nd=${d//_}
    grep -rl "$d" *.html | xargs sed -i "s/${d}/${nd}/g"
    mv $d $nd
done

cd ../..
cp -r _build/html/* ..
