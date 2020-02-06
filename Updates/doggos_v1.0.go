//Written by Dan Dochtermann
//Final Project: PCA and clustering of Dog breeds
//CMU Course 02-601, Fall 2018

package main

import (
	"fmt"
	"os"
	"bufio"
	"strings"
	"strconv"
)

// Struct containing all neccessary data for each individual dog within the dataset
type Pup struct {
	breed string 				//The breed specified within the dataset
	loci  map[string]float64 	//map of all STR loci of names to occurences
}

/* Misc notes/functioning code snippets

//runs after build with project.exe DogFiler.txt

	var input string
	fmt.Scanln(&input)
*/

func main () {
	//Checks for the correct number of command-line arguments
	if len(os.Args) != 2 {
		fmt.Println("Error: wrong number of input arguments")
		os.Exit(1)
	}
	
	inFilename := os.Args[1]
	pupMatrix := GeneratePupMatrix(inFilename)
}

func GeneratePupMatrix(inFile string) []Pup {
	var pupMatrix []Pup = make([]Pup, 0)
	var headers []string
	var locusNames []string
	fileLines := ReadFileLines(inFile)

	headers = strings.Split(fileLines[1], "	")
	locusNames = HeadersToLocusNames(headers)
	for i:=2; i<len(fileLines); i++ {
		pupMatrix = append(pupMatrix, GetNextPup(fileLines[i], locusNames))
	}

	return pupMatrix
}

// Opens the command-line input file and converts the data into a fileLines array of strings
// Homework 4 spatial.go was used as reference for this function format
func ReadFileLines(inFile string) []string{
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

func HeadersToLocusNames(headers []string) []string {
	var locusNames []string 
	for i:= 2; i<len(headers); i++ {
		if headers[i] != "" {
			locusNames = append(locusNames, headers[i])
		}
	}
	return locusNames
}

func GetNextPup(pupString string, locusNames []string) Pup {
	var newPup Pup 
	var pupInfo []string
	currentLoci := make(map[string]float64)

	pupInfo = strings.Split(pupString, "	")
	newPup.breed = pupInfo[1]

	for i:=2; i<len(pupInfo); i = i+2 {
		if pupInfo[i] != "" {
			var repeats float64
			a, erra := strconv.ParseFloat(pupInfo[i], 64)
			b, errb := strconv.ParseFloat(pupInfo[i+1], 64)
			if erra != nil || errb != nil {
				fmt.Println("Error converting STR repeats to float64")
				os.Exit(1)
			}
			repeats = (a+b)/2
			currentLoci[locusNames[(i/2)-1]] = repeats
		}

	}
	newPup.loci = currentLoci
	return newPup
}