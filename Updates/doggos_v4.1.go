//Written by Dan Dochtermann
//Final Project: PCA and clustering of Dog breeds
//CMU Course 02-601, Fall 2018

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

// Struct containing all neccessary data for each individual dog within the dataset
type Pup struct {
	breed string             //The breed specified within the dataset
	loci  map[string]float64 //map of all STR loci of names to occurences
	id    int                //id tag to recall later
}

// Will represent any 2D Matrix used
type DataMatrix [][]float64

type ClusterSet []Position

type TagInfo struct {
	label    []int
	position Position
}
type Position struct {
	x float64
	y float64
	z float64
}

/* Misc notes/functioning code snippets

//runs after build with [ project DogFiler.txt 2 ]
//Whippet, Longhaired is 2178
	var input string
	fmt.Scanln(&input)
*/

// Checks for the correct number of command-line arguments
// Returns effors when inputs are not recognized
// Calls functions to preform PCA and relevant algorithms
func main() {

	if len(os.Args) != 3 {
		fmt.Println("Error: wrong number of input arguments")
		fmt.Println("Input should be of the following format:")
		fmt.Println("[.exe name] [dataset filename] [PCA dimensions (2 or 3)]")
		os.Exit(1)
	}
	inFilename := os.Args[1]
	dimensions, err := strconv.Atoi(os.Args[2])
	if err != nil || dimensions > 3 {
		fmt.Println("Error parsing requested PCA dimensions, or too large a value entered")
		os.Exit(1)
	}

	maxClusterSize := 20
	//Get initial setup and breakdown once
	M, PupMatrix, locusNames := GeneratePupMatrix(inFilename)
	scores := PCA(M, dimensions)
	PupMatrices := MeanShiftAndCluster(PupMatrix, scores, dimensions, maxClusterSize)
	var FinalSet [][]Pup
	var tagArray []int
	var FinalTags [][]int
	RecursiveBreakdown(&FinalSet, &FinalTags, PupMatrices, locusNames, dimensions, maxClusterSize, tagArray)

	currentDepth := GetMaxDepth(FinalTags)
	var nextHsets []TagInfo
	var HSet []TagInfo

	for i := currentDepth; i > 0; i-- {
		HSet = append(GetSetsAtDepths(scores, FinalSet, FinalTags, i, dimensions), nextHsets...)
		nextHsets = AHC(HSet, i, dimensions)
	}
}

func GetMaxDepth(Tags [][]int) int {
	deepest := 0
	for i := range Tags {
		if len(Tags[i]) >= deepest {
			deepest = len(Tags[i])
		}
	}
	return deepest

}

// Agglomerative Heirarchal Clustering
// NEED TO RETURN MORE INFO, JUST COMPUNDING FOR NOW
func AHC(HSet []TagInfo, d int, dims int) []TagInfo {
	var toDo []int
	var x int
	var matches int
	var compoundedHSets []TagInfo
	//collect levels that will be grouped, //COMMETN ABOUT THIS PART HERE
	if d == 1 {
		compoundedHSets = append(compoundedHSets, HeirarchalCluster(HSet, d, dims))
	}
	if d != 1 {

		for i := range HSet {
			matches = 0
			x = HSet[i].label[d-2]
			for j := 0; j < len(toDo); j++ {
				if x == toDo[j] {
					matches++
				}
			}
			if matches == 0 {
				toDo = append(toDo, x)
			}
		}
		for i := range toDo {

			var toCluster []TagInfo
			for j := range HSet {
				if HSet[j].label[d-2] == toDo[i] {

					toCluster = append(toCluster, HSet[j])
				}
			}
			compoundedHSets = append(compoundedHSets, HeirarchalCluster(toCluster, d, dims))
		}
	}
	return compoundedHSets

}

//for now just make it compound them to make sure shit works
func HeirarchalCluster(HSet []TagInfo, d int, dims int) TagInfo {
	var newGroup TagInfo
	newTag := HSet[0].label[:d-1]
	var sumPos, newPos Position
	for i := range HSet {
		sumPos.x += HSet[i].position.x
		sumPos.y += HSet[i].position.y

		if dims == 3 {
			sumPos.z += HSet[i].position.z
		}
	}
	a := len(HSet)

	newPos.x = sumPos.x / float64(a)
	newPos.y = sumPos.y / float64(a)
	newPos.z = sumPos.z / float64(a)
	newGroup.position = newPos
	newGroup.label = newTag

	return newGroup
}

func GetSetsAtDepths(scores DataMatrix, pupSet [][]Pup, tagSet [][]int, currentDepth int, dims int) []TagInfo {
	var HSet []TagInfo
	for i := range tagSet {
		if len(tagSet[i]) == currentDepth {
			HSet = append(HSet, GetTagInfo(tagSet[i], scores, pupSet[i], dims))
		}
	}
	return HSet
}

func GetTagInfo(currentTag []int, scores DataMatrix, P []Pup, dims int) TagInfo {
	var NewSet TagInfo
	var newPos Position
	var sumPos Position
	for i := range P {
		j := P[i].id
		sumPos.x += scores[j][0]
		sumPos.y += scores[j][1]

		if dims == 3 {
			sumPos.z += scores[j][2]
		}
	}
	a := len(P)
	newPos.x = sumPos.x / float64(a)
	newPos.y = sumPos.y / float64(a)
	newPos.z = sumPos.z / float64(a)
	NewSet.position = newPos
	NewSet.label = currentTag

	return NewSet

}

// Recursively Breaks down clusters until they are smaller than the
// requested maxsize. Appends any clusters small enough to
// FinalSet, via pointers
func RecursiveBreakdown(FinalSet *[][]Pup, FinalTags *[][]int, PupMatrices [][]Pup, locusNames []string, dimensions int, maxSize int, tagArray []int) [][]Pup {
	for i := range PupMatrices {
		if len(PupMatrices[i]) <= maxSize {
			if len(PupMatrices[i]) > 0 {
				*FinalSet = append(*FinalSet, PupMatrices[i])
				sendArr := append(tagArray, i)
				*FinalTags = append(*FinalTags, sendArr)
				//fmt.Println(sendArr, "complete (", len(PupMatrices[i]), ")")
			}
		} else {
			tempArr := append(tagArray, i)
			//fmt.Println(tempArr, "being broken down (", len(PupMatrices[i]), ")")
			N := Create2DMatrix(PupMatrices[i], locusNames)
			scores := PCA(N, dimensions)
			NewPupMatrices := MeanShiftAndCluster(PupMatrices[i], scores, dimensions, maxSize)
			RecursiveBreakdown(FinalSet, FinalTags, NewPupMatrices, locusNames, dimensions, maxSize, tempArr)
		}
	}
	return *FinalSet
}

// Takes the command-line infile as an input
// Converts relevant filedata to an array of Pups
// Passes the pupMatrix on to a separate function
// Returns a DataMatrix of all relevant inFile data
func GeneratePupMatrix(inFile string) (DataMatrix, []Pup, []string) {
	var pupMatrix []Pup = make([]Pup, 0)
	fileLines := ReadFileLines(inFile)
	locusNames := HeadersToLocusNames(fileLines[1])

	for i := 2; i < len(fileLines); i++ {
		pupMatrix = append(pupMatrix, GetNextPup(fileLines[i], locusNames))
	}

	M := Create2DMatrix(pupMatrix, locusNames)

	return M, pupMatrix, locusNames
}

// Opens the command-line input file and converts the data into a fileLines array of strings
// Homework #4 spatial.go was used as reference for this function format
func ReadFileLines(inFile string) []string {
	datafile, err := os.Open(inFile)
	if err != nil {
		fmt.Println("error opening file")
		os.Exit(1)
	}

	fileLines := make([]string, 0)

	scanner := bufio.NewScanner(datafile)
	for scanner.Scan() {
		fileLines = append(fileLines, scanner.Text())
	}
	if scanner.Err() != nil {
		fmt.Println("error parsing file")
		os.Exit(1)
	}

	datafile.Close()
	return fileLines

}

// Takes the fileLine corresponding to header names
// splits them into individual headers, then processes those into an array of locus names
func HeadersToLocusNames(rawHeaders string) []string {
	var headers []string
	var locusNames []string
	headers = strings.Split(rawHeaders, "	")
	for i := 2; i < len(headers); i++ {
		if headers[i] != "" {
			locusNames = append(locusNames, headers[i])
		}
	}
	return locusNames
}

// takes a file line and the locusNames array
// splits the file line, assigning the breed to the newPup
// averages and converts the STR frequencies of both alleles
// into one value, and assigns the average to the map of that locus
// Returns the Pup with all data processed
func GetNextPup(pupString string, locusNames []string) Pup {
	var newPup Pup
	var pupInfo []string
	currentLoci := make(map[string]float64)

	pupInfo = strings.Split(pupString, "	")
	newPup.breed = pupInfo[1]
	id, err := strconv.Atoi(pupInfo[0])
	if err != nil {
		fmt.Println("Error parsing first column ID from dataset")
		os.Exit(1)
	}

	newPup.id = id - 1

	for i := 2; i < len(pupInfo); i = i + 2 {
		if pupInfo[i] != "" {
			var repeats float64
			a, erra := strconv.ParseFloat(pupInfo[i], 64)
			b, errb := strconv.ParseFloat(pupInfo[i+1], 64)
			if erra != nil || errb != nil {
				fmt.Println("Error converting STR repeats to float64")
				os.Exit(1)
			}
			repeats = (a + b) / 2
			currentLoci[locusNames[(i/2)-1]] = repeats
		}

	}
	newPup.loci = currentLoci
	return newPup
}

// converts a slice of Pups into a 2D matrix of standard format
// utilizes a locusNames array as it ranges across the PupMatrix
// in order to assure the same loci order is used for each Pup
func Create2DMatrix(p []Pup, locusNames []string) DataMatrix {
	M := make(DataMatrix, 0)
	for i := range p {
		n := make([]float64, 0)
		for j := range locusNames {
			n = append(n, p[i].loci[locusNames[j]])
		}
		M = append(M, n)
	}

	return M
}

// Takes a n by m DataMatrix as input
// returns it's tranformation, which sould be another
// DataMatrix of size m by n
func TransposeMatrix(M DataMatrix) DataMatrix {
	var T DataMatrix
	for i := 0; i < len(M[0]); i++ {
		newRow := make([]float64, len(M))
		T = append(T, newRow)
	}

	for i := 0; i < len(T); i++ {
		for j := 0; j < len(T[0]); j++ {
			T[i][j] = M[j][i]
		}
	}

	return T
}

// Takes two Matrices as input
// Performs standard matrix multiplication
// Returns the result of matrix multiplication
func MultiplyMatrix(M1, M2 DataMatrix) DataMatrix {

	var W DataMatrix
	for i := 0; i < len(M1); i++ {
		newRow := make([]float64, len(M2[0]))
		W = append(W, newRow)
	}
	for i := 0; i < len(W); i++ {
		for j := 0; j < len(W[i]); j++ {
			W[i][j] = DotProduct(M1, M2, i, j)
		}
	}

	return W
}

// Perfomrs the dot-product for MultiplyMatrix
// Takes two Matrices and the current row and column as input
// Creates arrays for the row of M1 and the column of M2
// Uses a third array to store their products
// Sums the third array and returns the value
func DotProduct(M1, M2 DataMatrix, i, j int) float64 {
	var a, b, c []float64
	for x := range M1[i] {
		a = append(a, M1[i][x])
	}
	for y := range M2 {
		b = append(b, M2[y][j])
	}
	for z := range a {
		c = append(c, (a[z] * b[z]))
	}

	prod := SumArray(c)

	return prod

}

// Sums contents across an array
// Takes a slice of float64 as input
// Returns the sum of it's contents
func SumArray(a []float64) float64 {
	var sum float64
	for i := range a {
		sum += a[i]
	}
	return sum
}

// Used for verifying results during testing
// Takes a DataMatrix as input
// Ranges across values and prints them to the console
func PrintMatrix(M DataMatrix) {
	for i := range M {
		for j := range M[i] {
			fmt.Print(M[i][j], " ")
		}
		fmt.Println()
	}
}

// Takes a DataMatrix and the requested dimensions as input
// Truncates to the first dims columns that will be used for PCA
// Returns the truncated DataMatrix
func TruncateMatrix(W DataMatrix, dims int) DataMatrix {
	var Wt DataMatrix
	for i := 0; i < len(W); i++ {
		var newRow []float64
		for j := 0; j < dims; j++ {
			newRow = append(newRow, W[i][j])
		}
		Wt = append(Wt, newRow)
	}
	return Wt
}

// Used for visualization of the PCA output
// Takes a DataMatrix and a given filename as input
// Prints them to a CSV file for plotting
func MatrixToCSV(filename string, M DataMatrix, P []Pup) {
	outfilename := filename + ".csv"
	outFile, err := os.Create(outfilename)
	if err != nil {
		fmt.Println("Error: couldn’t create", outfilename, "(is it currently open?)")
	}
	for i := 0; i < len(M); i++ {
		for j := 0; j < len(M[0]); j++ {
			fmt.Fprint(outFile, M[i][j], ",")
		}
		fmt.Fprint(outFile, P[i].breed)
		fmt.Fprintln(outFile)
	}
	outFile.Close()
}

func ClustersToCSV(filename string, C ClusterSet, dims int) {
	outfilename := filename + ".csv"
	outFile, err := os.Create(outfilename)
	if err != nil {
		fmt.Println("Error: couldn’t create", outfilename, "(is it currently open?)")
	}
	for i := 0; i < len(C); i++ {
		fmt.Fprint(outFile, C[i].x, ",", C[i].y, ",")
		if dims == 3 {
			fmt.Fprint(outFile, C[i].z, ",")
		}
		fmt.Fprintln(outFile)
	}
	outFile.Close()
}

// Takes a DataMatrix and the requested dimensions as input
// Calls all PCA relevant functions
// Returns a transformed and scaled DataMatrix
func PCA(PreM DataMatrix, dims int) DataMatrix {
	M := NormalizeMatrix(PreM)
	MT := TransposeMatrix(M)
	W := MultiplyMatrix(MT, M)
	Wt := TruncateMatrix(W, dims)
	scores := MultiplyMatrix(M, Wt)
	return scores
}

// Takes a DataMatrix as input
// Normalizes the data for each variable based on
// its average and standarad deviation
// Returns the normalized datamatrix
func NormalizeMatrix(M DataMatrix) DataMatrix {
	var MNorm DataMatrix
	for i := 0; i < len(M); i++ {
		newRow := make([]float64, len(M[0]))
		MNorm = append(MNorm, newRow)
	}
	for j := 0; j < len(M[0]); j++ {

		var a []float64
		var avg, stdev float64
		for i := 0; i < len(M); i++ {
			a = append(a, M[i][j])
		}
		avg, stdev = GetNormVars(a)
		for i := 0; i < len(M); i++ {
			MNorm[i][j] = Normalize(M[i][j], avg, stdev)
		}
	}
	return MNorm
}

// takes an array of a single variable's values as input
// calcualtes its average and standard deviation
// for normalization, returns avg and stdev
func GetNormVars(a []float64) (float64, float64) {
	var squares []float64
	numVals := len(a)
	avg := SumArray(a) / float64(numVals)
	for i := 0; i < len(a); i++ {
		squares = append(squares, (a[i]-avg)*(a[i]-avg))
	}
	stdev := math.Sqrt((SumArray(squares) / float64(numVals)))
	return avg, stdev
}

// takes a single value, the average and standard deviation
// normalizes it and returns the value
func Normalize(x float64, avg float64, stdev float64) float64 {
	return (x - avg) / stdev
}

// takes a DataMatrix, the input dimensions, and the maxClusterSize
// makes a field of initial clusters for MeanShift based on the
// maxSize and  the length of the input DataMatrix
// sets the cluster search radius equal to the largest dimension div
// returns the cluster field and the radius
func CreateClusters(M DataMatrix, dims int, maxSize int) (ClusterSet, float64) {
	var Clusters ClusterSet
	var zDiv, zStart float64
	d := (math.Sqrt(float64(len(M) / maxSize))) + 1
	minArr, maxArr := MinMax(M, dims)
	xDiv := (maxArr[0] - minArr[0]) / d
	xStart := minArr[0] + xDiv/float64(2)
	yDiv := (maxArr[1] - minArr[1]) / d
	yStart := minArr[1] + yDiv/float64(2)
	if dims == 3 {
		zDiv = (maxArr[2] - minArr[2]) / d
		zStart = minArr[0] + xDiv/float64(2)
	}
	for i := xStart; i < maxArr[0]; i += xDiv {
		for j := yStart; j < maxArr[1]; j += yDiv {
			if dims == 2 {
				var n Position
				n.x = i
				n.y = j
				Clusters = append(Clusters, n)
			}
			if dims == 3 {
				for k := zStart; k < maxArr[2]; k += zDiv {
					var n Position
					n.x = i
					n.y = j
					n.z = k
					Clusters = append(Clusters, n)
				}
			}
		}
	}
	var radius float64
	if xDiv >= yDiv {
		radius = xDiv
	} else {
		radius = yDiv
	}

	if dims == 3 && zDiv >= radius {
		radius = zDiv
	}
	return Clusters, radius
}

// Cluster method, shifts the clusters to areas of higher concetration
// takes the cluster field, current datamatrix, the radius and dims
// as input, moves all points within the cluster field over 100
// iterations, returns the an array of nulls where no data was within
// the radius, and the first cluster that wasn't empty
func (C ClusterSet) MeanShift(M DataMatrix, rad float64, dims int) ([]int, int) {
	var nulls []int
	firstNonNull := len(C) - 1
	for iterations := 0; iterations < 100; iterations++ {
		for a := range C {
			nearby := FindInRange(C[a], M, rad, dims)
			if len(nearby) != 0 {
				if a < firstNonNull {
					firstNonNull = a
				}
				avg := AvgPos(nearby, dims)
				C[a].x = avg.x
				C[a].y = avg.y
				C[a].z = avg.z
			} else {
				if iterations == 0 {
					nulls = append(nulls, a)
				}
			}
		}
	}
	return nulls, firstNonNull
}

// Takes a MeanShifted clusterfield and the first non-null, the radius
// the null array, and the dimensions
// Forms a new cluster field only consisting of non-duplicate values
// in addition to ignoring clusters that are within a half radius,
// essetially treating them as an equivalent cluster
// returns the new cluster field representing all means
func CleanUp(C ClusterSet, nulls []int, f int, rad float64, dims int) ClusterSet {
	var A ClusterSet
	var matches int
	var nullmatch bool
	var xDis2, yDis2, zDis2 float64
	A = append(A, C[f])
	for i := 0; i < len(C); i++ {
		matches = 0
		for j := 0; j < len(A); j++ {
			xDis2 = (C[i].x - A[j].x) * (C[i].x - A[j].x)
			yDis2 = (C[i].y - A[j].y) * (C[i].y - A[j].y)
			if dims == 3 {
				zDis2 = (C[i].z - A[j].z) * (C[i].z - A[j].z)
			}
			if math.Sqrt(xDis2+yDis2+zDis2) < rad/float64(2) {
				matches++
			}
		}
		if matches == 0 {
			nullmatch = false
			for k := 0; k < len(nulls); k++ {
				if nulls[k] == i {
					nullmatch = true
				}
			}
			if nullmatch == false {
				A = append(A, C[i])
			}
		}
	}
	return A
}

// finds the average position over an array of positions
func AvgPos(a []Position, dims int) Position {
	var xSum, ySum, zSum float64
	for i := range a {
		xSum += a[i].x
		ySum += a[i].y
		if dims == 3 {
			zSum += a[i].z
		}
	}
	var avg Position
	avg.x = xSum / float64(len(a))
	avg.y = ySum / float64(len(a))
	avg.z = zSum / float64(len(a))
	return avg
}

// searches the datamatrix for points that are within range of the
// moving cluster cetnter. If true, they are added to an array
// of positions, which is returned
func FindInRange(a Position, M DataMatrix, rad float64, dims int) []Position {
	var nearby []Position
	var xDis2, yDis2, zDis2 float64
	for i := 0; i < len(M); i++ {
		xDis2 = (a.x - M[i][0]) * (a.x - M[i][0])
		yDis2 = (a.y - M[i][1]) * (a.y - M[i][1])
		if dims == 3 {
			zDis2 = (a.z - M[i][2]) * (a.z - M[i][2])
		}
		if math.Sqrt(xDis2+yDis2+zDis2) <= rad {
			var newPos Position
			newPos.x = M[i][0]
			newPos.y = M[i][1]
			if dims == 3 {
				newPos.z = M[i][2]
			}
			nearby = append(nearby, newPos)
		}
	}
	return nearby

}

// Takes an integer x and returns 10^x
func Pow10(pow int) int {
	val := 1
	for i := 0; i < pow; i++ {
		val = val * 10
	}

	return val
}

// Finds the minimums and maximums in matrix
// returns 2 slices equal to the length of the input dimensions
// where each index is a max/min of that dimension (0=x, and so on)
func MinMax(M DataMatrix, dims int) ([]float64, []float64) {
	var mins = make([]float64, dims)
	var maxs = make([]float64, dims)

	for i := 0; i < dims; i++ {
		var tempMax, tempMin float64
		tempMax = -999999999
		tempMin = 999999999
		for j := 0; j < len(M); j++ {
			if M[j][i] >= tempMax {
				tempMax = M[j][i]
			}
			if M[j][i] <= tempMin {
				tempMin = M[j][i]
			}
		}
		mins[i] = tempMin
		maxs[i] = tempMax
	}
	return mins, maxs
}

// Calls all related functions to run Mean shift clustering
// Creates a clusterfield, meanshifts them, cleans them up
// and splits them into a new set of matrices
// when that set is of length 1, meaning it only made one
// cluster, it is run a second time with the detection
// radius cut in half so clusters can be distinct
// returns the new PupMatrices

// [May want to implement heirarchal from here, then more from branches]
// Either that or somehow maintain the level at which each breakdown
// happens
func MeanShiftAndCluster(M []Pup, scores DataMatrix, dims int, maxSize int) [][]Pup {
	var PupMatrices [][]Pup

	C, rad := CreateClusters(scores, dims, maxSize)
	nulls, firstNonNull := C.MeanShift(scores, rad, dims)
	cores := CleanUp(C, nulls, firstNonNull, rad, dims)
	PupMatrices = SplitToMatrices(M, cores, scores, dims)

	if len(PupMatrices) == 1 {
		C, rad := CreateClusters(scores, dims, maxSize)
		rad = rad / float64(2)
		nulls, firstNonNull := C.MeanShift(scores, rad, dims)
		cores := CleanUp(C, nulls, firstNonNull, rad, dims)
		PupMatrices = SplitToMatrices(M, cores, scores, dims)
	}

	return PupMatrices
}

// Does the splitting part of MeanShiftAndCluster
// Ranges over the input slice, finding the closest cluster
// in the input set of cores. Copies that Pup from
// the input to the index of the output slice corresponding
// to the core it is closest to
func SplitToMatrices(M []Pup, cores ClusterSet, scores DataMatrix, dims int) [][]Pup {
	var PupMatrices = make([][]Pup, len(cores))
	//var scoreSplit = make([]DataMatrix, len(cores))
	var j int
	for i := range M {
		var newPup Pup
		//get closest cluster index in cores
		j = ClosestCluster(cores, scores, dims, i)
		newPup = M[i]
		PupMatrices[j] = append(PupMatrices[j], newPup)
		//scoreSplit[j] = append(scoreSplit[j], scores[i])
	}

	return PupMatrices
}

// takes i, the index of original input matrix
// calculates the distance to all of the available cluster cores
// returns nearest core
func ClosestCluster(cores ClusterSet, scores DataMatrix, dims int, i int) int {
	var xDis2, yDis2, zDis2 float64
	var distance float64
	distance = 999999999
	var nearest int
	for j := 0; j < len(cores); j++ {
		xDis2 = (cores[j].x - scores[i][0]) * (cores[j].x - scores[i][0])
		yDis2 = (cores[j].y - scores[i][1]) * (cores[j].y - scores[i][1])
		if dims == 3 {
			zDis2 = (cores[j].z - scores[i][2]) * (cores[j].z - scores[i][2])
		}
		if math.Sqrt(xDis2+yDis2+zDis2) <= distance {
			distance = math.Sqrt(xDis2 + yDis2 + zDis2)
			nearest = j
		}

	}
	return nearest
}
