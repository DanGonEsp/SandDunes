--[[
-------------------------------
Save this file as simplecsv.lua
-------------------------------

Example use: Suppose file csv1.txt is:

1.23,70,hello
there,9.81,102
x,y,,z
8,1.243,test

Then the following
-------------------------------------
local csvfile = require "simplecsv"
local m = csvfile.read('./csv1.txt') -- read file csv1.txt to matrix m
print(m[2][3])                       -- display element in row 2 column 3 (102)
m[1][3] = 'changed'                  -- change element in row 1 column 3
m[2][3] = 123.45                     -- change element in row 2 column 3
csvfile.write('./csv2.txt', m)       -- write matrix to file csv2.txt
-------------------------------------

will produce file csv2.txt with the contents:

1.23,70,changed
there,9.81,123.45
x,y,,z
8,1.243,test

the read method takes 4 parameters:
path: the path of the CSV file to read - mandatory
sep: the separator character of the fields. Optionsl, defaults to ','
tonum: whether to convert fields to numbers if possible. Optional. Defaults to true
null: what value should null fields get. Optional. defaults to ''
]]

module(..., package.seeall)

---------------------------------------------------------------------
function string:split(sSeparator, nMax, bRegexp)
    if sSeparator == '' then
        sSeparator = ','
    end

    if nMax and nMax < 1 then
        nMax = nil
    end

    local aRecord = {}

    if self:len() > 0 then
        local bPlain = not bRegexp
        nMax = nMax or -1

        local nField, nStart = 1, 1
        local nFirst,nLast = self:find(sSeparator, nStart, bPlain)
        while nFirst and nMax ~= 0 do
            aRecord[nField] = self:sub(nStart, nFirst-1)
            nField = nField+1
            nStart = nLast+1
            nFirst,nLast = self:find(sSeparator, nStart, bPlain)
            nMax = nMax-1
        end
        aRecord[nField] = self:sub(nStart)
    end

    return aRecord
end

---------------------------------------------------------------------
function read(path, sep, tonum, null)
    tonum = tonum or true
    sep = sep or ','
    null = null or ''
    local csvFile = {}
    local file = assert(io.open(path, "r"))
    for line in file:lines() do
        fields = line:split(sep)
        if tonum then -- convert numeric fields to numbers
            for i=1,#fields do
                local field = fields[i]
                if field == '' then
                    field = null
                end
                fields[i] = tonumber(field) or field
            end
        end
        table.insert(csvFile, fields)
    end
    file:close()
    
    local num_rows = #csvFile					-- Output: Number of rows:
	local num_cols = #csvFile[1]				-- Output: Number of columns:
	print("num_rows = "..num_rows .. "  num_cols = " ..num_cols)
	
	for i_ind = 2, num_rows do
		for j_ind = 1, num_cols do
			if (type(csvFile[i_ind][j_ind]) == nil or type(csvFile[i_ind][j_ind]) ~= "number") then
				print ("Invalid Input parametrs (ERROR) value[" .. i_ind-1 .."][".. j_ind .."]  = " .. type(csvFile[i_ind][j_ind])); exit();
			end
		end
	end

    return csvFile,num_rows,num_cols
    
    
end

---------------------------------------------------------------------
function write(path, data, sep)
    sep = sep or ','
    local file = assert(io.open(path, "w"))
    for i=1,#data do
        for j=1,#data[i] do
            if j>1 then file:write(sep) end
            file:write(data[i][j])
        end
        file:write('\n')
    end
    file:close()
end





---------------------------------------------------------------------
function MergeDebugPVD(path, name)
    -- Helper to prepend 'cd "path" &&' to each command
    local function run(cmd)
        local full_cmd = string.format('cd "%s" && %s', path, cmd)
        local ok = os.execute(full_cmd)
        if ok ~= 0 then
            print("Command failed:", full_cmd)
        end
    end

    -- Commands translated from your shell script
    -- 1) Update timestep in .pvtu files using Perl
    run([[for f in *.pvtu; do perl -i -0777 -pe 'if(/_iter(\d+)/){$i=$1; s/(<Time timestep=")\d+(")/$1$i$2/}' "$f"; done]])
	print("------------- 10%")
    -- 2) Rename .pvtu files
    run([[for f in *.pvtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvtu/\1_\3_iter\2.pvtu/')"; done]])
	print("------------- 20%")
    -- 3) Rename .pvd files
    run([[for f in *.pvd; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvd/\1_\3_iter\2.pvd/')"; done]])
	print("------------- 30%")
    -- 4) Rename .vtu files
    run([[for f in *.vtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.vtu/\1_\3_iter\2.vtu/')"; done]])
	print("------------- 40%")
    -- 5) Update file paths inside .pvd files
    run([[for f in *.pvd; do perl -i -pe 's/(<DataSet .*file=".*?)(?:_iter(\d+))(_.*?)(\.pvtu")/$1$3_iter$2$4/' "$f"; done]])
	print("------------- 50%")
    -- 6) Update timestep in .pvd files
    run([[for f in *.pvd; do perl -i -pe 's/(timestep=")\d+(".*file="[^"]*_iter(\d+)[^"]*")/$1.($3+0).$2/e' "$f"; done]])
	print("------------- 60%")
    -- 7) Update Piece sources in .pvtu files
    run([[for f in *.pvtu; do perl -i -pe 's/(<Piece Source=".*?)(_iter\d+)(.*)(\.vtu")/$1$3$2$4/' "$f"; done]])
end


--[[
	os.execute("
	cd path
	for f in path/*.pvtu; do perl -i -0777 -pe 'if(/_iter(\d+)/){$i=$1; s/(<Time timestep=")\d+(")/$1$i$2/}' "$f"; done
	for f in *.pvtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvtu/\1_\3_iter\2.pvtu/')"; done
	for f in *.pvd; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvd/\1_\3_iter\2.pvd/')"; done
	for f in *.vtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.vtu/\1_\3_iter\2.vtu/')"; done
	for f in *.pvd; do perl -i -pe 's/(<DataSet .*file=".*?)(?:_iter(\d+))(_.*?)(\.pvtu")/$1$3_iter$2$4/' "$f"; done
	for f in *.pvd; do perl -i -pe 's/(timestep=")\d+(".*file="[^"]*_iter(\d+)[^"]*")/$1.($3+0).$2/e' "$f"; done
	for f in *.pvtu; do perl -i -pe 's/(<Piece Source=".*?)(_iter\d+)(.*)(\.vtu")/$1$3$2$4/' "$f"; done
	")
    
    

	for f in *.pvtu; do perl -i -0777 -pe 'if(/_iter(\d+)/){$i=$1; s/(<Time timestep=")\d+(")/$1$i$2/}' "$f"; done
	for f in *.pvtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvtu/\1_\3_iter\2.pvtu/')"; done
	for f in *.pvd; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvd/\1_\3_iter\2.pvd/')"; done
	for f in *.vtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.vtu/\1_\3_iter\2.vtu/')"; done
	for f in *.pvd; do perl -i -pe 's/(<DataSet .*file=".*?)(?:_iter(\d+))(_.*?)(\.pvtu")/$1$3_iter$2$4/' "$f"; done
	for f in *.pvd; do perl -i -pe 's/(timestep=")\d+(".*file="[^"]*_iter(\d+)[^"]*")/$1.($3+0).$2/e' "$f"; done
	for f in *.pvtu; do perl -i -pe 's/(<Piece Source=".*?)(_iter\d+)(.*)(\.vtu")/$1$3$2$4/' "$f"; done


python3 - <<'EOF'
import xml.etree.ElementTree as ET, glob, re

def merge_pvd(pattern, output_name):
    master = ET.Element("VTKFile", type="Collection", version="0.1")
    col = ET.SubElement(master, "Collection")

    def get_iter(f):
        m = re.search(r'_iter(\d+)', f)
        return int(m.group(1)) if m else 0

    lines = []
    for fname in glob.glob(pattern):
        try:
            tree = ET.parse(fname)
            root = tree.getroot()
            for ds in root.findall(".//DataSet"):
                lines.append((get_iter(ds.get("file")), ET.tostring(ds, encoding="unicode")))
        except ET.ParseError:
            print(f"Skipping invalid XML: {fname}")

    for _, ds_xml in sorted(lines, key=lambda x: x[0]):
        col.append(ET.fromstring(ds_xml))

    ET.ElementTree(master).write(output_name, encoding="utf-8", xml_declaration=True)
    print(f"Created {output_name} with pattern {pattern}")



	# Merge Defect files
	merge_pvd("*_Defect_*.pvd", "master_Defect.pvd")

	# Merge Solution files
	merge_pvd("*_Solution_*.pvd", "master_Solution.pvd")

	# Merge Correction files
	merge_pvd("*_Correction_*.pvd", "master_Correction.pvd")

EOF

    

    

# Merge Defect files
merge_pvd("*_Defect_*.pvd", "master_Defect_"..name..".pvd")

# Merge Solution files
merge_pvd("*_Solution_*.pvd", "master_Solution_"..name..".pvd")

# Merge Correction files
merge_pvd("*_Correction_*.pvd", "master_Correction_"..name..".pvd")
    
end
]]














--[[


	for f in *.pvtu; do perl -i -0777 -pe 'if(/_iter(\d+)/){$i=$1; s/(<Time timestep=")\d+(")/$1$i$2/}' "$f"; done
	for f in *.pvtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvtu/\1_\3_iter\2.pvtu/')"; done
	for f in *.pvd; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.pvd/\1_\3_iter\2.pvd/')"; done
	for f in *.vtu; do mv "$f" "$(echo "$f" | sed -E 's/(.*)_iter([0-9]{3})_(.*)\.vtu/\1_\3_iter\2.vtu/')"; done
	for f in *.pvd; do perl -i -pe 's/(<DataSet .*file=".*?)(?:_iter(\d+))(_.*?)(\.pvtu")/$1$3_iter$2$4/' "$f"; done
	for f in *.pvd; do perl -i -pe 's/(timestep=")\d+(".*file="[^"]*_iter(\d+)[^"]*")/$1.($3+0).$2/e' "$f"; done
	for f in *.pvtu; do perl -i -pe 's/(<Piece Source=".*?)(_iter\d+)(.*)(\.vtu")/$1$3$2$4/' "$f"; done


python3 - <<'EOF'
import xml.etree.ElementTree as ET, glob, re

def merge_pvd(pattern, output_name):
    master = ET.Element("VTKFile", type="Collection", version="0.1")
    col = ET.SubElement(master, "Collection")

    def get_iter(f):
        m = re.search(r'_iter(\d+)', f)
        return int(m.group(1)) if m else 0

    lines = []
    for fname in glob.glob(pattern):
        try:
            tree = ET.parse(fname)
            root = tree.getroot()
            for ds in root.findall(".//DataSet"):
                lines.append((get_iter(ds.get("file")), ET.tostring(ds, encoding="unicode")))
        except ET.ParseError:
            print(f"Skipping invalid XML: {fname}")

    for _, ds_xml in sorted(lines, key=lambda x: x[0]):
        col.append(ET.fromstring(ds_xml))

    ET.ElementTree(master).write(output_name, encoding="utf-8", xml_declaration=True)
    print(f"Created {output_name} with pattern {pattern}")



	# Merge Defect files
	merge_pvd("*_Defect_*.pvd", "master_Defect.pvd")

	# Merge Solution files
	merge_pvd("*_Solution_*.pvd", "master_Solution.pvd")

	# Merge Correction files
	merge_pvd("*_Correction_*.pvd", "master_Correction.pvd")




EOF


]]



