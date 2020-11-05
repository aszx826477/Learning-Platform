from pyspark import SparkContext

logFile = "README.md"  # Should be some file on your system
sc = SparkContext()
logData = sc.textFile(logFile)

numAs = logData.filter(lambda s: 'a' in s).count()
numBs = logData.filter(lambda s: 'b' in s).count()

print ("Lines with a:", numAs, " lines with b: ", numBs)

sc.stop()
