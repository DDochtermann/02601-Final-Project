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
}

// Will represent any 2D Matrix used
type DataMatrix [][]float64

type ClusterSet []Position

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

	M := GeneratePupMatrix(inFilename)
	scores := PCA(M, dimensions)
	C, rad := CreateClusters(scores, dimensions)
	A := MeanShift(C, scores, rad, dimensions)
	ClustersToCSV("Cleanedupclusters", A, dimensions)
}

// Takes the command-line infile as an input
// Converts relevant filedata to an array of Pups
// Passes the pupMatrix on to a separate function
// Returns a DataMatrix of all relevant inFile data
func GeneratePupMatrix(inFile string) DataMatrix {
	var pupMatrix []Pup = make([]Pup, 0)
	fileLines := ReadFileLines(inFile)
	locusNames := HeadersToLocusNames(fileLines[1])

	for i := 2; i < len(fileLines); i++ {
		pupMatrix = append(pupMatrix, GetNextPup(fileLines[i], locusNames))
	}

	M := Create2DMatrix(pupMatrix, locusNames)

	return M
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
func MatrixToCSV(filename string, M DataMatrix) {
	outfilename := filename + ".csv"
	outFile, err := os.Create(outfilename)
	if err != nil {
		fmt.Println("Error: couldn’t create", outfilename, "(is it currently open?)")
	}
	for i := 0; i < len(M); i++ {
		for j := 0; j < len(M[0]); j++ {
			fmt.Fprint(outFile, M[i][j], ",")
		}
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

func CreateClusters(M DataMatrix, dims int) (ClusterSet, float64) {
	var Clusters ClusterSet
	var zDiv, zStart float64
	minArr, maxArr := MinMax(M, dims)
	xDiv := (maxArr[0] - minArr[0]) / float64(10)
	xStart := minArr[0] + xDiv/float64(2)
	yDiv := (maxArr[1] - minArr[1]) / float64(10)
	yStart := minArr[1] + yDiv/float64(2)
	if dims == 3 {
		zDiv = (maxArr[2] - minArr[2]) / float64(10)
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

func MeanShift(C ClusterSet, M DataMatrix, rad float64, dims int)  ClusterSet {
	var nulls []int
	for iterations := 0; iterations < 100; iterations++ {
		for a := range(C) {
			nearby := FindInRange(C[a], M, rad, dims)
			if len(nearby) != 0 {
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
	A := CleanUp(C, nulls)
	return A

}

func CleanUp(C ClusterSet, nulls []int) ClusterSet {
	var A ClusterSet
	var matches int
	var nullmatch bool
	A = append(A, C[0])
	for i:=0; i<len(C); i++ {
		matches = 0
		for j:=0; j<len(A); j++ {
			if C[i].x == A[j].x && C[i].y == A[j].y && C[i].z == A[j].z {
				matches++
			}
		}
		if matches == 0 {
			nullmatch = false
			for k := 0; k<len(nulls); k++ {
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

func AvgPos(a []Position, dims int) Position {
	var xSum, ySum, zSum float64
	for i:=range(a) {
		xSum += a[i].x
		ySum += a[i].y
		if dims == 3 {
			zSum += a[i].z
		}
	}
	var avg Position
	avg.x = xSum/float64(len(a))
	avg.y = ySum/float64(len(a))
	avg.z = zSum/float64(len(a))
	return avg
}

func FindInRange(a Position, M DataMatrix, rad float64, dims int) []Position {
	var nearby []Position
	var xDis2, yDis2, zDis2 float64
	for i:=0; i<len(M); i++ {
		xDis2 = (a.x - M[i][0])*(a.x - M[i][0])
		yDis2 = (a.y - M[i][1])*(a.y - M[i][1])
		if dims == 3 {
			zDis2 = (a.z - M[i][2])*(a.z - M[i][2])
		}
		if math.Sqrt(xDis2 + yDis2 + zDis2) <= rad {
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

func Pow10(pow int) int {
	val := 1
	for i := 0; i < pow; i++ {
		val = val * 10
	}

	return val
}


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
