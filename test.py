import time


import cnstools_2.file_handlers as fhs
import source.cnstools.filetype_classes as fcs

new_t = 0
for i in range(10):
    t = time.time()
    new = fhs.gff3.Handler("test_files/Gmax.gff")
    new.to_bed("test_files/new.gff",type_list=["CDS"],genome="Gmax")
    new_t+=(time.time()-t)
print new_t/10.0

old_t = 0
for i in range(10):
    t = time.time()
    old = fcs.Gff3("test_files/Gmax.gff")
    old.to_bed(type_list=["CDS"],genome="Gmax").save_file("test_files/old.gff")
    old_t+=(time.time()-t)
print old_t/10.0

