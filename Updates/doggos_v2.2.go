//Written by Dan Dochtermann
//Final Project: PCA and clustering of Dog breeds
//CMU Course 02-601, Fall 2018

package main

import (
	"bufio"
	"fmt"
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
type DataMatrix [][] float64

/* Misc notes/functioning code snippets

//runs after build with project DogFiler.txt
//Whippet, Longhaired is 2178

	var input string
	fmt.Scanln(&input)
*/

func main() {
	//Checks for the correct number of command-line arguments
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
	
	MatrixToCSV("scores", scores)
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
	for i := 0; i<len(M[0]); i++ {
		newRow := make([]float64, len(M))
		T = append(T, newRow)
	} 

	for i := 0; i<len(T); i++ {
		for j:=0; j<len(T[0]); j++ {
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
	for i := 0; i<len(M1); i++ {
		newRow := make([]float64, len(M2[0]))
		W = append(W, newRow)
	}
	for i:=0; i<len(W); i++ {
		for j:=0; j<len(W[i]); j++ {
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
func DotProduct(M1, M2 DataMatrix, i, j int) float64{
	var a, b, c []float64
	for x := range M1[i] {
		a = append(a, M1[i][x])
	}
	for y := range M2 {
		b = append(b, M2[y][j])
	}
	for z := range a {
		c = append(c, (a[z]*b[z]))
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
		for j:= range M[i] {
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
	for i:=0; i<len(W); i++ {
		var newRow []float64
		for j:=0; j<dims; j++ {
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
	for i:=0; i<len(M); i++ {
		for j:= 0; j<len(M[0]); j++ {
			fmt.Fprint(outFile, M[i][j], ",")
		}
		fmt.Fprintln(outFile)
	}
	outFile.Close()
}

func PCA(M DataMatrix, dims int) DataMatrix {
	MT := TransposeMatrix(M)
	W := MultiplyMatrix(MT, M)
	Wt := TruncateMatrix(W, dims)
	scores := MultiplyMatrix(M, Wt)
	return scores
}