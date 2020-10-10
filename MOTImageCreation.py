
import os, numpy, imageio, random, copy, math
from scipy.spatial import distance

########################## MOT Square Image Creation ###############################
# This function creates images with various squares moving in random directions for MOT conditions.
# Input parameters for MOTSquaresImageCreation:
    # numImages - number of images to be generated corresponding to different movements of squares
    # numSquares - number of sqaures in the image
    # imageRows - number of rows of pixels in the image
    # imageColumns - number of columns of pixels in the image
    # squareSize - squareSize+1 will be the size in number of pixels of a square
    # squareMoveStepSize - How many pixels will the center of the square will move each time
    
def MOT_Squares_ImageCreation(numImages, numSquares, imageRows, imageColumns, squareSize, squareMoveStepSize):
    dirPath = os.getcwd() + "/Stimuli/" + str("MOT_Squares") + "/"
    outImageStimulus = []
    squareCentersX = []
    squareCentersY = []
    availablePositionsCentral = []
    availablePositionsCentral2 = []
    #random.seed(1000)
    # Not including corner (squareSize/2) positions in availablePositions list
    for i in range(int(squareSize / 2), int(imageRows - squareSize / 2)):
        for j in range(int(squareSize / 2), int(imageColumns - squareSize / 2)):
            availablePositionsCentral.append(int(i * imageColumns + j))

    for i in range(int(squareSize / 2 + squareMoveStepSize), int(imageRows - squareSize / 2 - squareMoveStepSize)):
        for j in range(int(squareSize / 2 + squareMoveStepSize), int(imageColumns - squareSize / 2 - squareMoveStepSize)):
            availablePositionsCentral2.append(int(i * imageColumns + j))

    centerXPrevious = []  # list of x positions of the centers of all the squares in the previous move
    centerYPrevious = []  # list of y positions of the centers of all the squares in the previous move

    for h in range(numImages):
        if h == 0:
            centerX = [0] * numSquares
            centerY = [0] * numSquares
            im = numpy.zeros(shape=(imageRows, imageColumns))
            availablePositions = copy.deepcopy(availablePositionsCentral)
            for s in range(numSquares):
                center = random.randint(0, len(availablePositions))
                centerX[s] = availablePositions[center] // imageColumns
                centerY[s] = availablePositions[center] % imageColumns

                for i in range(int(-squareSize / 2), int(squareSize / 2 + 1)):
                    for j in range(int(-squareSize / 2), int(squareSize / 2 + 1)):
                        im[centerX[s] + i][centerY[s] + j] = 255

                for i in range(-squareSize, squareSize + 1):
                    for j in range(-squareSize, squareSize + 1):
                        q = (centerX[s] + i) * imageColumns + (centerY[s] + j)
                        if q in availablePositions:
                            availablePositions.remove(q)

            centerXPrevious = centerX
            centerYPrevious = centerY
            squareCentersX.append(centerX)
            squareCentersY.append(centerY)

        else:
            centerX = [0] * numSquares
            centerY = [0] * numSquares
            im = numpy.zeros(shape=(imageRows, imageColumns))
            availablePositions = copy.deepcopy(availablePositionsCentral2)
            for s in range(numSquares):
                availableMoveDirections = list(range(8))
                currentPosition = 0
                while currentPosition not in availablePositions:
                    #print("length = " + str(len(availableMoveDirections)-1))
                    index = random.randint(0, (len(availableMoveDirections) - 1))
                    #print("index = "+str(index))
                    moveDirection = availableMoveDirections[index]
                    if moveDirection == 0:
                        centerX[s] = centerXPrevious[s] + squareMoveStepSize
                        centerY[s] = centerYPrevious[s]
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 1:
                        centerX[s] = centerXPrevious[s] + squareMoveStepSize
                        centerY[s] = centerYPrevious[s] + squareMoveStepSize
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 2:
                        centerX[s] = centerXPrevious[s]
                        centerY[s] = centerYPrevious[s] + squareMoveStepSize
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 3:
                        centerX[s] = centerXPrevious[s] - squareMoveStepSize
                        centerY[s] = centerYPrevious[s] + squareMoveStepSize
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 4:
                        centerX[s] = centerXPrevious[s] - squareMoveStepSize
                        centerY[s] = centerYPrevious[s]
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 5:
                        centerX[s] = centerXPrevious[s] - squareMoveStepSize
                        centerY[s] = centerYPrevious[s] - squareMoveStepSize
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 6:
                        centerX[s] = centerXPrevious[s]
                        centerY[s] = centerYPrevious[s] - squareMoveStepSize
                        availableMoveDirections.remove(moveDirection)
                    elif moveDirection == 7:
                        centerX[s] = centerXPrevious[s] + squareMoveStepSize
                        centerY[s] = centerYPrevious[s] - squareMoveStepSize
                        availableMoveDirections.remove(moveDirection)

                    currentPosition = int(centerX[s] * imageColumns + centerY[s])

                for i in range(int(-squareSize / 2), int(squareSize / 2 + 1)):
                    for j in range(int(-squareSize / 2), int(squareSize / 2 + 1)):
                        im[centerX[s] + i][centerY[s] + j] = 255

                for i in range(-squareSize, squareSize + 1):
                    for j in range(-squareSize, squareSize + 1):
                        q = (centerX[s] + i) * imageColumns + (centerY[s] + j)
                        if q in availablePositions:
                            availablePositions.remove(q)

            centerXPrevious = centerX
            centerYPrevious = centerY
            squareCentersX.append(centerX)
            squareCentersY.append(centerY)
        
        im = im.astype('uint8')
        filename = dirPath + "MOT_Squares_" + str(h) + ".bmp"
        if not os.path.exists(dirPath):
            os.makedirs(dirPath)
        imageio.imwrite(filename, im)
        outImageStimulus.append(im)

    #print(squareCentersX)
    file = dirPath + "MOTSquares Moving Squares.GIF"
    imageio.mimwrite(file, outImageStimulus, duration=1.0)
    return squareCentersX, squareCentersY

# Running the MOT square image creation code
'''
numImages = 10
imageRows = 86
imageColumns = 86
numSquares = 5
squareSize = 10
squareMoveStepSize = 5
MOT_Squares_ImageCreation(numImages, numSquares, imageRows, imageColumns, squareSize, squareMoveStepSize)
'''

######################### MOT PS1988 #########################
# This simulates experiemnt 1 from Pylyshyn and Storm 1988
# 10 plus signs with 1 to 5 targets, never collided with each other,
# velocity and direction (among 8) changed every few 100 ms, reflected from edges
# Input parameters for MOTSquaresImageCreation:
	# numImages - number of images to be generated corresponding to different movements of pluses
	# numCross - number of plus signs in the image
	# imageRows - number of rows of pixels in the image
	# imageColumns - number of columns of pixels in the image
	# crossSize - will be the size in number of pixels in one direction of plus sign

def MOT_PS1988_ImageCreation(totalNumImages, numCross, imageRows, imageColumns, crossSize, minDistance, minSpeed, maxSpeed, fixationSquareSize, fps, secondsPerFrame, endProbeFlashingNoFrames):
    dirPath = os.getcwd() + "/Stimuli/" + str("MOT_PS1988") + "/"
    print("totalnumImages = "+str(totalNumImages))
    moveNumImages = totalNumImages - endProbeFlashingNoFrames # Number of frames for which the targets move
    print("numMovingImages = " + str(moveNumImages))

    outImageStimulus = []
    crossCentersX = []
    crossCentersY = []
    availablePositionsCentral = []

    # Not including corner (crossSize/2) positions in available positions
    for i in range(int((crossSize-1) / 2+1), int(imageRows - (crossSize-1) / 2)):
        for j in range(int((crossSize-1) / 2+1), int(imageColumns - (crossSize-1) / 2)):
            availablePositionsCentral.append(int(i * imageColumns + j))

    fixation = [] # contains a list of all position values of the fixation square
    for i in range(int((imageRows + 1) / 2 - (fixationSquareSize - 1) / 2), int((imageRows + 1) / 2 + (fixationSquareSize - 1) / 2 + 1)):
        for j in range(int((imageRows + 1) / 2 - (fixationSquareSize - 1) / 2), int((imageRows + 1) / 2 + (fixationSquareSize - 1) / 2 + 1)):
            q = int(i * imageColumns + j)
            fixation.append(q)
    availablePositionsCentral = [x for x in availablePositionsCentral if x not in fixation]# Removing the center fixation square positions from available positions

    #centerXPrevious = []  # list of x positions of the centers of all the squares in the previous move
    #centerYPrevious = []  # list of y positions of the centers of all the squares in the previous move
    futureCentersX = []
    futureCentersY = []
    h = 0
    initialOverlap = 1
    while h < (moveNumImages-1): # h th image has already been created in the previous loop except for when h=0
        changeFrame = int(random.choice([(fps * 0.1 * x) for x in range(1,6)])) # Number of frames after which to change direction and velocity of crosses, calculated using fps and changes every few 100 miliseconds
        print("h = " + str(h))
        if h == 0:# Initializing the crosses positions randomly
            availablePositions = copy.deepcopy(availablePositionsCentral)
            while initialOverlap == 1:
                centerXtemp = []
                centerYtemp = []
                for s in range(numCross):
                    center = random.randint(0, len(availablePositions)) # Randomly choosing a position for one cross
                    centerXtemp.append(availablePositions[center] // imageColumns) # Converting the position value to X coordinate
                    centerYtemp.append(availablePositions[center] % imageColumns) # Converting the position value to Y coordinate

                # Checking if the image's crosses are within minDistance from each other
                coordinates = [[i, j] for (i, j) in zip(centerXtemp, centerYtemp)]
                dist = distance.cdist(coordinates, coordinates, 'euclidean')
                dist[dist == 0] = numpy.amax(dist)  # Replacing the zeroes i.e. the values of distance from themselves to the maximum value in the array
                if numpy.any(dist < minDistance):
                    initialOverlap = 1
                    print("Breaking from too close crosses initially")
                else:
                    initialOverlap = 0
                    centerXPrevious = centerXtemp # Converting the position value to X coordinate
                    centerYPrevious = centerYtemp # Converting the position value to Y coordinate
                    futureCentersX.append(centerXPrevious) # to create the initial image, adding the first image to the future images list
                    futureCentersY.append(centerYPrevious) # to create the initial image, adding the first image to the future images list
                    crossCentersX.append(centerXPrevious)
                    crossCentersY.append(centerYPrevious)

        # Creating centers for next changeFrame number of iterations and making sure the crosses dont come close to each other
        overlap = 1
        counter = 0
        while overlap == 1:
            # So that the if the code gets stuck we can restart from 1st image itself
            counter = counter + 1
            if counter >= 5000:
                print("Breaking counter")
                h = 0
                break

            direction = [random.randint(0, 7) for _ in range(numCross)]  # Determine new direction for each cross, will be a list of size equal to number of crosses
            velocity = [random.randint(minSpeed, maxSpeed) for _ in range(numCross)]  # Determine new velocity for each cross, will be a list of size equal to number of crosses
            trialX = [] # List of list of number of changeFrame amount of frames with numCross centers
            trialY = [] # List of list of number of changeFrame amount of frames with numCross centers
            centerXPreviousLocal = copy.deepcopy(centerXPrevious) # So we keep centerXPrevious in memory so we can go back to this if overlap happens
            centerYPreviousLocal = copy.deepcopy(centerYPrevious)

            f = 0
            while f < changeFrame:
                # Find the X and Y coordinates of centers for new image's crosses
                x = []
                y = []
                crossesTogether = []  # Contains a single list of crosses positions
                s = 0
                while s < numCross:
                    if direction[s] == 0:
                        crossX = int(round(centerXPreviousLocal[s] + velocity[s] * secondsPerFrame))
                        crossY = int(round(centerYPreviousLocal[s]))
                    elif direction[s] == 1:
                        crossX = int(round(centerXPreviousLocal[s] + velocity[s] * secondsPerFrame/math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] + velocity[s] * secondsPerFrame/math.sqrt(2)))
                    elif direction[s] == 2:
                        crossX = int(round(centerXPreviousLocal[s]))
                        crossY = int(round(centerYPreviousLocal[s] + velocity[s] * secondsPerFrame))
                    elif direction[s] == 3:
                        crossX = int(round(centerXPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] + velocity[s] * secondsPerFrame / math.sqrt(2)))
                    elif direction[s] == 4:
                        crossX = int(round(centerXPreviousLocal[s] - velocity[s] * secondsPerFrame))
                        crossY = int(round(centerYPreviousLocal[s]))
                    elif direction[s] == 5:
                        crossX = int(round(centerXPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))
                    elif direction[s] == 6:
                        crossX = int(round(centerXPreviousLocal[s]))
                        crossY = int(round(centerYPreviousLocal[s] - velocity[s] * secondsPerFrame))
                    elif direction[s] == 7:
                        crossX = int(round(centerXPreviousLocal[s] + velocity[s] * secondsPerFrame / math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))

                    # Saving coordinates of crosses
                    oneCross = []
                    oneCrossX = []
                    oneCrossY = []
                    for i in range(int(-(crossSize - 1) / 2), int((crossSize - 1) / 2 + 1)):
                        oneCross.append(int((crossX + i) * imageColumns + crossY))
                        oneCrossX.append(int(crossX + i))
                        oneCrossY.append(int(crossY))
                    for j in range(int(-(crossSize - 1) / 2), int((crossSize - 1) / 2 + 1)):
                        oneCross.append(int(crossX * imageColumns + crossY + j))
                        oneCrossX.append(int(crossX))
                        oneCrossY.append(int(crossY + j))
                    #oneCross = list(set(oneCross)) # Removing reccurring pixel indices

                    # If any cross are going outside the background we are bouncing them from edges i.e. changing their direction
                    #pixelsIndex = list(range(0, imageRows * imageColumns))
                    #outsideCrossPixel = [i for i in oneCross if i not in pixelsIndex]  # finding the crosses which are outside the background
                    outside = [i for i in oneCrossX if i >= imageRows or i < 0]
                    outside.extend([i for i in oneCrossY if i >= imageColumns or i < 0])
                    #print("len(outsideCrossPixel) = " + str(len(outsideCrossPixel)))
                    if len(outside) > 0:
                        #print("Outside Pixels = "+str(outsideCrossPixel))
                        # Change direction of these crosses, bounce off edges
                        if direction[s] <= 3:
                            direction[s] = direction[s] + 4
                        elif direction[s] >= 4:
                            direction[s] = direction[s] - 4
                        continue

                    crossesTogether.extend(oneCross)
                    x.append(crossX)
                    y.append(crossY)
                    s = s + 1

                # Checking if any cross has overlapped with the center fixation square
                if any(f in crossesTogether for f in fixation):
                    overlap = 1
                    #print("Breaking from fixation overlap")
                    break

                # Checking if the new image's crosses are within minDistance from each other
                # if overlap happens (overlap = 1) - break for loop
                # if overlap = 0 for this value of f - replace centerXPreviousLocal = x to create the next image in line
                coordinates = [[i,j] for (i,j) in zip(x,y)]
                dist = distance.cdist(coordinates, coordinates, 'euclidean')
                dist[dist == 0] = numpy.amax(dist) # Replacing the zeroes i.e. the values of distance from themselves to the maximum value in the array
                if numpy.any(dist < minDistance):
                    overlap = 1
                    print("Breaking from too close crosses")
                    break
                else:
                    overlap = 0

                centerXPreviousLocal = copy.deepcopy(x)
                centerYPreviousLocal = copy.deepcopy(y)

                trialX.append(x)
                trialY.append(y)
                f = f + 1

        futureCentersX.extend(trialX)
        futureCentersY.extend(trialY)
        centerXPrevious = futureCentersX[len(futureCentersX) - 1] # Save the previous crosses center's positions for next image creation
        centerYPrevious = futureCentersY[len(futureCentersY) - 1]
        h = h + changeFrame  # Update the number of images whose crosses center have been determined

    print("Finding centers done\n Creating images...")
    ### Creating all the images
    outputCentersX = []
    outputCentersY = []
    for locNum in range(len(futureCentersX)):
        im = numpy.zeros(shape=(imageRows, imageColumns))  # Empty image array

        # Creating the center white fixation square
        for q in fixation:
            fixX = q // imageColumns
            fixY = q % imageColumns
            im[fixX][fixY] = 255

        # Creating the crosses
        for s in range(numCross):
            for i in range(int(-(crossSize-1) / 2), int((crossSize-1) / 2 + 1)):
                im[futureCentersX[locNum][s] + i][futureCentersY[locNum][s]] = 255
            for j in range(int(-(crossSize-1) / 2), int((crossSize-1) / 2 + 1)):
                im[futureCentersX[locNum][s]][futureCentersY[locNum][s] + j] = 255

        # Creating .bmp of each frame of stimulus
        im = im.astype('uint8')
        #print("creating image number = " + str(locNum))
        filename = dirPath + "MOT_PS1988_" + str(locNum) + ".bmp"
        if not os.path.exists(dirPath):
            os.makedirs(dirPath)
        if locNum < moveNumImages: # The image creating stops when the number of images created are numImages
            imageio.imwrite(filename, im)
            outImageStimulus.append(im) # Saving in a list to create gif of moving crosses
            outputCentersX.append(futureCentersX[locNum])
            outputCentersY.append(futureCentersY[locNum])

    # Adding extra frames with last images
    for t in range(endProbeFlashingNoFrames):
        filename = dirPath + "MOT_PS1988_" + str(moveNumImages + t) + ".bmp"
        im = outImageStimulus[-1]
        imageio.imwrite(filename, im)
        outImageStimulus.append(im)
        outputCentersX.append(outputCentersX[-1])
        outputCentersY.append(outputCentersY[-1])

    #print("outputCentersX = "+str(outputCentersX))

    print("number of Images Created = "+str(len(outImageStimulus)))
    # Creating gif of moving crosses
    print("done\nCreating gif of moving crosses")
    file = dirPath + "MOT PS1988 Moving Crosses.GIF"
    imageio.mimwrite(file, outImageStimulus, duration=secondsPerFrame)
    return outputCentersX, outputCentersY


# Running the PS1988 image creation code
experimentRunTime = 7*1000#random.randint(7,15)
fps = 30  # Frames per second, should be a multiple of 10 so that frames changees at constant time intervals
secondsPerFrame = 1.0 / fps  # in seconds
endProbeFlashingNoFrames = 60
numImages = int(fps * experimentRunTime / 1000) + endProbeFlashingNoFrames
numObjects = 10
numTargets = 5
imageRows = 113  # in pixels, needs to be odd to add the center fixation square
imageColumns = 113  # in pixels
objectSize = 7 #13 # in pixels, 1 center pixel and half size on each 4 sides of the center
fixationSquareSize = 5 #13 # in pixels, size of center concetration square
minDistance = 11 #37 # in pixels, minimum distance between objects center to center - 24(end-end minimum) + 13(size of object)
minSpeed = 30 #40 # in pixels/sec, Minimum speed of object movement
maxSpeed = 49
X, Y = MOT_PS1988_ImageCreation(numImages, numObjects, imageRows,
                                imageColumns, objectSize, minDistance,
                                minSpeed, maxSpeed,
                         fixationSquareSize, fps, secondsPerFrame,
                                endProbeFlashingNoFrames)


######################### MOT Yantis1992 #########################
# This simulates experiemnt 1 from Yantis 1992
# 10 plus signs with 3 to 5 targets, never collided with each other,
# velocity and direction (among 8) changed every trajectory duration choosen in the beginning for each cross, reflected from edges
# Input parameters for MOTSquaresImageCreation:
	# numImages - number of images to be generated corresponding to different movements of pluses
	# numCross - number of plus signs in the image
	# imageRows - number of rows of pixels in the image
	# imageColumns - number of columns of pixels in the image
	# crossSize - will be the size in number of pixels in one direction of plus sign
    # cushion - number of pixels from the center of the cross or fixation square within which no other object can enter
    # speed - speed of cross movement in pixels/second
    # fixationSquareSize - size in pixels of the central fixation square
    # fps - frames per second
    # secondsPerFrame - seconds one frame is on

def MOT_Yantis1992_ImageCreation(numImages, numCross, imageRows, imageColumns, crossSize, cushion, speed, fps, secondsPerFrame):
    dirPath = os.getcwd() + "/Stimuli/" + str("MOT_Yantis1992") + "/"
    outImageStimulus = []
    crossCentersX = []
    crossCentersY = []
    availablePositionsCentral = []
    #print("numImages = "+str(numImages))
    #print("seconds = "+str(secondsPerFrame))
    # Not including corner (crossSize/2) positions in available positions
    for i in range(int((crossSize-1) / 2+1), int(imageRows - (crossSize-1) / 2)):
        for j in range(int((crossSize-1) / 2+1), int(imageColumns - (crossSize-1) / 2)):
            availablePositionsCentral.append(int(i * imageColumns + j))

    fixation = [] # contains positions of the fixation square and minDistance and half of crossSize from fixation square center
    fixDist = cushion + (crossSize - 1) / 2 # Distance from center which contains fixation square and its cushion and crosses cushion
    for i in range(int((imageRows + 1) / 2 - fixDist), int((imageRows + 1) / 2 + fixDist + 1)):
        for j in range(int((imageColumns + 1) / 2 - fixDist), int((imageColumns + 1) / 2 + fixDist + 1)):
            q = int(i * imageColumns + j)
            fixation.append(q)
    availablePositionsCentral = [x for x in availablePositionsCentral if x not in fixation]# Removing the center fixation square positions from available positions

    centerXPrevious = []  # list of x positions of the centers of all the squares in the previous move
    centerYPrevious = []  # list of y positions of the centers of all the squares in the previous move
    futureCentersX = []
    futureCentersY = []
    # Remove later
    trajectoryDuration = int(random.choice(range(210, 840, 30)))*speed  # Number of frames after which to change direction and velocity of crosses
    h = 0
    while h < (numImages-1): # h th image has already been created in the previous loop except for when h=0
        changeFrame = int(random.choice([(fps * 0.1 * x) for x in range(1,6)])) # Number of frames after which to change direction and velocity of crosses, calculated using fps and changes every few 100 miliseconds
        print("h = " + str(h))
        #print("changeFrame = "+str(changeFrame))
        if h == 0:# Initializing the crosses positions randomly
            availablePositions = copy.deepcopy(availablePositionsCentral)
            for s in range(numCross):
                center = random.randint(0, len(availablePositions)) # Randomly choosing a position for one cross
                centerXPrevious.append(availablePositions[center] // imageColumns) # Converting the position value to X coordinate
                centerYPrevious.append(availablePositions[center] % imageColumns) # Converting the position value to Y coordinate
            futureCentersX.append(centerXPrevious) # to create the initial image, adding the first image to the future images list
            futureCentersY.append(centerYPrevious) # to create the initial image, adding the first image to the future images list
            crossCentersX.append(centerXPrevious)
            crossCentersY.append(centerYPrevious)

        # Creating centers for next changeFrame number of iterations and making sure the crosses dont come close to each other
        overlap = 1
        while overlap == 1:
            direction = [random.randint(0, 7) for _ in range(numCross)]  # Determine new direction for each cross, will be a list of size equal to number of crosses
            velocity = [random.randint(minSpeed, maxSpeed) for _ in range(numCross)]  # Determine new velocity for each cross, will be a list of size equal to number of crosses
            #print("Direction = "+str(direction))
            #print("velocity = "+str(velocity))
            trialX = [] # List of list of number of changeFrame amount of frames with numCross centers
            trialY = [] # List of list of number of changeFrame amount of frames with numCross centers
            centerXPreviousLocal = copy.deepcopy(centerXPrevious) # So we keep centerXPrevious in memory so we can go back to this if overlap happens
            centerYPreviousLocal = copy.deepcopy(centerYPrevious)

            f = 0
            while f < changeFrame:
                #print ("F = " + str(f))
            #for f in range(changeFrame):
                # Find the X and Y coordinates of centers for new image's crosses
                x = []
                y = []
                crossesTogether = []  # Contains a single list of crosses positions
                s = 0
                while s < numCross:
                    #print("s = " + str(s))
                    #print("centerXPreviousLocal = " + str(centerXPreviousLocal))
                    if direction[s] == 0:
                        crossX = int(round(centerXPreviousLocal[s] + velocity[s] * secondsPerFrame))
                        crossY = int(round(centerYPreviousLocal[s]))
                    elif direction[s] == 1:
                        crossX = int(round(centerXPreviousLocal[s] + velocity[s] * secondsPerFrame/math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] + velocity[s] * secondsPerFrame/math.sqrt(2)))
                    elif direction[s] == 2:
                        crossX = int(round(centerXPreviousLocal[s]))
                        crossY = int(round(centerYPreviousLocal[s] + velocity[s] * secondsPerFrame))
                    elif direction[s] == 3:
                        crossX = int(round(centerXPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] + velocity[s] * secondsPerFrame / math.sqrt(2)))
                    elif direction[s] == 4:
                        crossX = int(round(centerXPreviousLocal[s] - velocity[s] * secondsPerFrame))
                        crossY = int(round(centerYPreviousLocal[s]))
                    elif direction[s] == 5:
                        crossX = int(round(centerXPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))
                    elif direction[s] == 6:
                        crossX = int(round(centerXPreviousLocal[s]))
                        crossY = int(round(centerYPreviousLocal[s] - velocity[s] * secondsPerFrame))
                    elif direction[s] == 7:
                        crossX = int(round(centerXPreviousLocal[s] + velocity[s] * secondsPerFrame / math.sqrt(2)))
                        crossY = int(round(centerYPreviousLocal[s] - velocity[s] * secondsPerFrame / math.sqrt(2)))

                    # Saving coordinates of crosses
                    oneCross = []
                    oneCrossX = []
                    oneCrossY = []
                    for i in range(int(-(crossSize - 1) / 2), int((crossSize - 1) / 2 + 1)):
                        oneCross.append(int((crossX + i) * imageColumns + crossY))
                        oneCrossX.append(int(crossX + i))
                        oneCrossY.append(int(crossY))
                    for j in range(int(-(crossSize - 1) / 2), int((crossSize - 1) / 2 + 1)):
                        oneCross.append(int(crossX * imageColumns + crossY + j))
                        oneCrossX.append(int(crossX))
                        oneCrossY.append(int(crossY + j))
                    #oneCross = list(set(oneCross)) # Removing reccurring pixel indices

                    # If any cross are going outside the background we are bouncing them from edges i.e. changing their direction
                    #pixelsIndex = list(range(0, imageRows * imageColumns))
                    #outsideCrossPixel = [i for i in oneCross if i not in pixelsIndex]  # finding the crosses which are outside the background
                    outside = [i for i in oneCrossX if i >= imageRows or i < 0]
                    outside.extend([i for i in oneCrossY if i >= imageColumns or i < 0])
                    #print("len(outsideCrossPixel) = " + str(len(outsideCrossPixel)))
                    if len(outside) > 0:
                        #print("Outside Pixels = "+str(outsideCrossPixel))
                        # Change direction of these crosses, bounce off edges
                        if direction[s] <= 3:
                            direction[s] = direction[s] + 4
                        elif direction[s] >= 4:
                            direction[s] = direction[s] - 4
                        continue

                    crossesTogether.extend(oneCross)
                    x.append(crossX)
                    y.append(crossY)
                    s = s + 1

                # Checking if any cross has overlapped with the center fixation square
                if any(f in crossesTogether for f in fixation):
                    overlap = 1
                    #print("Breaking from fixation overlap")
                    break

                # Checking if the new image's crosses are within the cushion of each other
                # if overlap happens (overlap = 1) - break for loop
                # if overlap = 0 for this value of f - replace centerXPreviousLocal = x to create the next image in line
                coordinates = [[i,j] for (i,j) in zip(x,y)]
                dist = distance.cdist(coordinates, coordinates, 'euclidean')
                dist[dist == 0] = numpy.amax(dist) # Replacing the zeroes i.e. the values of distance from themselves to the maximum value in the array
                if numpy.any(dist < cushion):
                    overlap = 1
                    #print("Breaking from too close crosses")
                    break
                else:
                    overlap = 0

                centerXPreviousLocal = copy.deepcopy(x)
                centerYPreviousLocal = copy.deepcopy(y)

                trialX.append(x)
                trialY.append(y)
                f = f + 1

        futureCentersX.extend(trialX)
        futureCentersY.extend(trialY)
        centerXPrevious = futureCentersX[len(futureCentersX) - 1] # Save the previous crosses center's positions for next image creation
        centerYPrevious = futureCentersY[len(futureCentersY) - 1]
        h = h + changeFrame  # Update the number of images whose crosses center have been determined

    #print("FutureCentersX = " + str(futureCentersX))
    #print("FutureCentersY = " + str(futureCentersY))
    print("Finding centers done")
    ### Creating all the images
    #print("len(futureCentersX) = " + str(len(futureCentersX)))
    for locNum in range(len(futureCentersX)):
        im = numpy.zeros(shape=(imageRows, imageColumns))  # Empty image array

        # Creating the center white fixation square
        for q in fixation:
            fixX = q // imageColumns
            fixY = q % imageColumns
            im[fixX][fixY] = 255

        # Creating the crosses
        for s in range(numCross):
            for i in range(int(-(crossSize-1) / 2), int((crossSize-1) / 2 + 1)):
                im[futureCentersX[locNum][s] + i][futureCentersY[locNum][s]] = 255
            for j in range(int(-(crossSize-1) / 2), int((crossSize-1) / 2 + 1)):
                im[futureCentersX[locNum][s]][futureCentersY[locNum][s] + j] = 255

        # Creating .bmp of each frame of stimulus
        im = im.astype('uint8')
        #print("creating image number = " + str(locNum))
        filename = dirPath + "MOT_Yantis1992_" + str(locNum) + ".bmp"
        if not os.path.exists(dirPath):
            os.makedirs(dirPath)
        if locNum < numImages: # The image creating stops when the number of images created are numImages
            imageio.imwrite(filename, im)
            outImageStimulus.append(im) # Saving in a list to create gif of moving crosses

    # Creating gif of moving crosses
    print("Creating gif of moving crosses")
    file = dirPath + "MOT Yantis1992 Moving Crosses.GIF"
    imageio.mimwrite(file, outImageStimulus, duration=secondsPerFrame)
    return futureCentersX, futureCentersY

'''
# Running the Yantis1992 image creation code
experimentRunTime = 7.5 * 1000  # in milliseconds
milliSecondsPerFrame = 30  # in milliseconds
secondsPerFrame = float(milliSecondsPerFrame) / 1000
fps = round(1.0/secondsPerFrame, 1)  # Frames per second, should be a multiple of 10 so that frames changes at constant time intervals
numImages = 250
numCross = 10
imageRows = 117  # in pixels, needs to be odd to add the center fixation square
imageColumns = 83  # in pixels
crossSize = 7 # in pixels, 1 center pixel and half size on each 4 sides of the center
fixationSquareSize = 5  # in pixels, size of center concetration square
minDistance = 14  # in pixels, minimum distance between objects center to center - 24(end-end minimum) + 13(size of object)
speed = 49  # in pixels/sec, Maximum  speed of object movement

X, Y = MOT_Yantis1992_ImageCreation(numImages, numCross, imageRows, imageColumns, crossSize, minDistance, speed,
                         fixationSquareSize, fps, secondsPerFrame)
print("X = "+str(X))
print("Y = "+str(Y))
'''