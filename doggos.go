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

// Global Variable, the string that will be complied for the Newick output
var HeirarchalString string

// Struct containing all neccessary data for each individual dog within the dataset
type Pup struct {
	breed string             // The breed specified within the dataset
	loci  map[string]float64 // map of all STR loci of names to occurences
	id    int                // id tag to recall later
}

// Trees and Nodes to be used for Phylogenetic tree construction
type Tree []*Node

type Node struct {
	label    string  // The label of the node for consistent numbering 
					 // and verification of Mean Shift visualization
	children []*Node // The collection of child nodes in an array
	score 	 Position// To be used for construction of a heirarchal tree
	weight 	 int 	 // also for heirarchal tree, will contain the amount
					 // of data points for calculation
	branch	float64  // represents the branch length to the node
	id 		int  	 // used for verification in test outputs

}

// Will represent any 2D Matrix used for PCA purposes
type DataMatrix [][]float64

// To represent a field of moving centroids in mean-shift clustering
type ClusterSet []Position

// For convenience of tracking 2 or 3D positions
type Position struct {
	x float64
	y float64
	z float64
}

// Checks for the correct number of command-line arguments
// Returns errors when inputs are not recognized
// Calls functions to preform PCA and relevant algorithms
func main() {
	if len(os.Args) != 2 {
		fmt.Println("Error: wrong number of input arguments")
		fmt.Println("Input should be of the following format:")
		fmt.Println("[.exe name] [dataset filename]")
		os.Exit(1)
	}

	// Default values for normal run
	inFilename := os.Args[1]
	dimensions := 3
	maxClusterSize := 20
	newickOut := "Newick"

	// Allows for modified PCA dimensions and clustering values for testing
	if inFilename == "testmode" {
		var testfile string
		fmt.Println("testmode detected.")
		fmt.Println("Please enter the input filename:")
		fmt.Scanln(&testfile)

		var testdims string
		fmt.Println("Please enter the requested PCA dimensions (2 or 3):")
		fmt.Scanln(&testdims)
		intdims, err := strconv.Atoi(testdims)
		if err != nil || intdims > 3 || intdims < 2 {
			fmt.Println("Error parsing requested PCA dimensions, invalid number entered")
			os.Exit(1)
		}

		var testsize string
		fmt.Println("Please enter the requested max cluster size:")
		fmt.Scanln(&testsize)
		intsize, err := strconv.Atoi(testsize)
		if err != nil {
			fmt.Println("Error parsing requested PCA dimensions, invalid number entered. Exiting.")
			os.Exit(1)
		}

		var testOut string
		fmt.Println("Please enter the Newick output filename (without the file extension):")
		fmt.Scanln(&testOut)

		inFilename = testfile
		dimensions = intdims
		maxClusterSize = intsize
		newickOut = testOut
	}

	// Get initial setup and breakdown once
	M, PupMatrix, locusNames := GeneratePupMatrix(inFilename)
	// Run Principal Component Analysis and MS cluster once
	scores := PCA(M, dimensions)
	PupMatrices := MeanShiftAndCluster(PupMatrix, scores, dimensions, maxClusterSize)
	
	// Will eventually hold the final clusters
	var FinalSet [][]Pup

	// Breaks down clusters over the maxClusterSize
	RecursiveBreakdown(&FinalSet, PupMatrices, locusNames, dimensions, maxClusterSize)
	
	// Creates the Heirarchal Cluster Tree
	HCT := CreateHCTree(FinalSet, scores, dimensions)

	// Prints the resulting tree to a text file in Newick format
	HeirarchalString += "("
	PrintHCTreeToNewick(HCT, HCT[len(HCT)-1])
	HeirarchalString += ");"
	NewickToTXT(newickOut)

}

// Takes the command-line infile as an input
// Converts relevant filedata to an array of Pups
// Passes the pupMatrix on to a separate function
// Returns a DataMatrix of all relevant inFile data
func GeneratePupMatrix(inFile string) (DataMatrix, []Pup, []string) {
	var pupMatrix []Pup = make([]Pup, 0)
	fileLines := ReadFileLines(inFile)
	locusNames := HeadersToLocusNames(fileLines[1])

	// Starting from the third line of the file, collect all samples
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
		fmt.Println("Error opening file, please verify filename and directory")
		os.Exit(1)
	}

	// Make an array of strings, add each line found to it
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
	// Marker names are every two columns, collect them
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

	// Split each line by detecting tabs in between values
	pupInfo = strings.Split(pupString, "	")
	newPup.breed = pupInfo[1]
	id, err := strconv.Atoi(pupInfo[0])
	if err != nil {
		fmt.Println("Error parsing first column ID from dataset")
		os.Exit(1)
	}

	// Offset ID by 1 to match other 0 indexed calls
	newPup.id = id - 1

	// Average the two STR values for each marker into one
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
// (Unused in final build, but important for testing)
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
// (Unused in final build, but important for testing)
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
	// Get largest radius
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
	var j int
	for i := range M {
		var newPup Pup
		//get closest cluster index in cores
		j = ClosestCluster(cores, scores, dims, i)
		newPup = M[i]
		PupMatrices[j] = append(PupMatrices[j], newPup)
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

// Recursively Breaks down clusters until they are smaller than the
// requested maxsize. Appends any clusters small enough to
// FinalSet, via pointers
func RecursiveBreakdown(FinalSet *[][]Pup, PupMatrices [][]Pup, locusNames []string, dimensions int, maxSize int) [][]Pup {
	for i := range PupMatrices {
		// Exit condition, if the set is under max size
		if len(PupMatrices[i]) <= maxSize {
			if len(PupMatrices[i]) > 0 {
				*FinalSet = append(*FinalSet, PupMatrices[i])
			}
		} else {
			// Recursion, break down sets larger than max size
			N := Create2DMatrix(PupMatrices[i], locusNames)
			scores := PCA(N, dimensions)
			NewPupMatrices := MeanShiftAndCluster(PupMatrices[i], scores, dimensions, maxSize)
			RecursiveBreakdown(FinalSet, NewPupMatrices, locusNames, dimensions, maxSize)
		}
	}
	return *FinalSet
}

// Works nearly the same as MatrixToCSV, but instead takes
// a clusterset as input, printing it to a file for
// plot visualization and testing purposes
// (Unused in final build, but important for testing)
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

// Takes the final Tree and the origin node as original input
// Recursively loops through when child nodes are not empty
// Fills the Newick string according to the standard format
func PrintHCTreeToNewick(T Tree, parent *Node) {
	if len(parent.children) != 0 {
		HeirarchalString += "("
	}
	for i:= range(parent.children) {
		PrintHCTreeToNewick(T, parent.children[i])
		if i == 0 {
			HeirarchalString += ","
		}
		if i == 1 {
			HeirarchalString += ")"
		}
	}
	HeirarchalString += parent.label
	HeirarchalString += ":"
	HeirarchalString += fmt.Sprintf("%F", parent.branch)
}


// Creates and opens an output file
// Stores the NewickString to it for upload to
// phlogeny displaying websites
func NewickToTXT(filename string) {
	outfilename := filename + ".txt"
	outFile, err := os.Create(outfilename)
	if err != nil {
		fmt.Println("Error: couldn’t create", outfilename, "(is it currently open?)")
	}
	fmt.Fprint(outFile, HeirarchalString)
	outFile.Close()
	fmt.Println("Output generated to", outfilename)
}

// Modified from UPGMA
// Initializes the HC tree, assigning labels, scores and weights
// To terminal nodes, leaves internal nodes unlabeled/unscored/unweighted
func CreateHCTree(PupSets [][]Pup, scores DataMatrix, dims int) Tree {
	numPups := len(PupSets)
	var t Tree = make([]*Node, 2*numPups-1)
	for i := range t {
		var vx Node
		if i < numPups { 
			vx.label = MakeStringForHCLabel(PupSets[i])
			vx.weight = len(PupSets[i])
			vx.score = AvgPos(DMtoPos(PupSets[i], scores, dims), dims)
			vx.id = i
		}
		t[i] = &vx
	}
	t.GroupNodes(numPups, dims)
	return t
}

// Ranges over the Pups in one cluster
// Makes a map and prints it to a string so the final output format
// Is easy to understand (Pug x1 /Chihuahua x3, etc.)
func MakeStringForHCLabel(P []Pup) string {
	m := make(map[string]int)
	for j := range P {
		m[P[j].breed]++
	}
	var n string
	for i := range m {
		n += i + " x" + strconv.Itoa(m[i]) + " /"
	}
	return n
}

// Takes a datamatrix and the Pups that correspond to it and retrieves
// their coordinates in the datamatrix, appending them to a new slice
// of x/y/z positions
func DMtoPos(P []Pup, scores DataMatrix, dims int) []Position {
	var a []Position
	for i:=range(P) {
		var newPos Position
		newPos.x = scores[P[i].id][0]
		newPos.y = scores[P[i].id][1]
		if dims == 3 {
			newPos.z = scores[P[i].id][2]
		}
		a = append(a, newPos)
	}
	return a
}

// Tree method, groups nodes based on heirarchal clustering
// finds the two closest clusters, creates a new node, 
// assigns that nodes values, and points the next available
// node in the tree to it.
// makes use of toGroup array to keep track of which nodes
// are available for grouping
func (t Tree) GroupNodes(numLeaves int, dims int) {
	var toGroup []int
	var a, b int
	var branch float64
	for i:=0; i<numLeaves; i++ {
		toGroup = append(toGroup, i)
	}

	for len(toGroup) > 1 {
		a, b, branch = FindClosestClusters(t, toGroup, dims)
		var newNode Node
		newNode.children = []*Node{t[toGroup[a]], t[toGroup[b]]}
		s, w := CalcScoreAndWeight(t[toGroup[a]],t[toGroup[b]])
		newNode.score = s
		newNode.weight = w
		newNode.id = numLeaves
		// Set newNode to the next empty node (same value as numLeaves,
		// which is now our counter variable)
		t[numLeaves] = &newNode
		t[toGroup[a]].branch = branch/2
		t[toGroup[b]].branch = branch/2 
		toGroup = append(toGroup, numLeaves)
		numLeaves++
		// move the last two indexes of the toGroup array to the 
		// positions of a and b, then crop off the last two indeces
		// effectively deleting a and b from the array
		toGroup[a] = toGroup[len(toGroup)-1]
		toGroup[b] = toGroup[len(toGroup)-2]
		toGroup = toGroup[:len(toGroup)-2]
	}

}

// searches within the tree among the nodes available for grouping, and 
// returns the two closest clusters as well as their distance to one 
// another for branch values
func FindClosestClusters(HCT Tree, toGroup []int, dims int) (int, int, float64) {
	var minDist float64
	var a, b int
	minDist = 999999999999999
	for i:=0; i<len(toGroup); i++ {
		for j:=i+1; j<len(toGroup); j++ {
			var d float64
			d = CalcDistance(HCT[toGroup[i]].score, HCT[toGroup[j]].score, dims)
			if d <= minDist {
				minDist = d
				a = i
				b = j
			}
		}
	}
	return a, b, minDist
}

// Simply calculates the distance between two 3D points
func CalcDistance(a Position, b Position, dims int) float64 {
	var xDis2, yDis2, zDis2 float64
	xDis2 = (a.x-b.x)*(a.x-b.x)
	yDis2 = (a.y-b.y)*(a.y-b.y)
	if dims == 3 {
		zDis2 = (a.z-b.z)*(a.z-b.z)
	}
	return math.Sqrt(xDis2+yDis2+zDis2)
}

// Takes two nodes and calculates a weighted average position for them
// returns that position as well as the sum of their weights for the parent node
func CalcScoreAndWeight(a, b *Node) (Position, int) {
	newWeight := a.weight + b.weight
	newX := (a.score.x * float64(a.weight) + b.score.x * float64(b.weight))/float64(newWeight)
	newY := (a.score.y * float64(a.weight) + b.score.y * float64(b.weight))/float64(newWeight)
	newZ := (a.score.z * float64(a.weight) + b.score.z * float64(b.weight))/float64(newWeight)
	var Pos Position
	Pos.x, Pos.y, Pos.z = newX, newY, newZ

	return Pos, newWeight
}

// The end. =]