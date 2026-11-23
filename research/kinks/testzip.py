import sys
import zipfile

z = zipfile.ZipFile(sys.argv[1])
for info in z.infolist():
    print(info.filename, info.file_size, info.compress_size)
print("testzip:", z.testzip())
