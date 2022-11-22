%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Farid Javadnejad (javadnef@oregonstate.edu)
% Created: 02/19/2016
% Modified: 06/01/2017
%
% FLIR E6

% This program gets the Interior and Exterior parameters of RGB cameras
% used for generating a SfM 3D model and the point cloud itsled. Then
% add a thermal field to the point cloud by integrating the data from
% thermal images acquired siminatously with oiptical images. A projection
% model is required to define the relatioship between RGB and thermal
% cameras.
%
% The required files are an XML file that included IO parametes (image
% size, focal length, principal point, distortion coefffincetns, etc.) for
% the RGB camera.
% The EO paramets are the SfM solution for location and orientetaiton of
% each camera in a coordainte system called "block coordiante system (CS)".
% Block coordainte system is an arbitrary CS defined intternally in SfM
% which does not have metric content. The XML file defines the
% transfoamtion paramewters from transforminng form chcunk CS to a real
% world coordainte.
%
% The EO and IO papremeters are used to iterate though an image-based
% reconstructed  3D point cloud and estimate the image coordinates [u v]
% for the projection of the 3D point with [Xw Yw Zw] in the images that
% are reporjected from
%
% Later, the pixel value for matching point on the pair thermal images in
% extarcted and is added to the ponit cloud as a new field.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initiate
clear all;
format long g;

% TO_DISPLAY = 1;
TO_APPLY_MASK = 1; %to apply a mask on thermal images

% Thermal Readings while the scale bar was fixed
TEMP_MAX_READING = 31.6; %celsius
TEMP_MIN_READING = 18.5; %celsius

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILES AND FOLDERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('X:\Farid\Publications\PERS Paper\MATLAB\');

% folders
WORKING_FOLDER = strcat('Y:\Farid\Projects\2017-07-11 Kearney Thermal\Thermal Processing\');
FOLDER_RGB_IMAGES = 'Y:\Farid\Projects\2017-07-11 Kearney Thermal\Images\RGB\';
FOLDER_TIR_IMAGES = 'Y:\Farid\Projects\2017-07-11 Kearney Thermal\Images\TIR\';

% pointcloud files reading
fileCloudRead = 'Y:\Farid\Projects\2017-07-11 Kearney Thermal\SfM Results\kearney_dense.txt'; % X,Y,Z,R,G,B
% output pointcloud file
fileCloudWrite = [fileCloudRead(1:end-4) '_TIR_3D_myRT.txt'];

% camera files
fileCameraRGB_XML = [WORKING_FOLDER 'kearney_dense.xml']; % IO, EO, Tw, Rw
fileThermalCameraCalib = [WORKING_FOLDER 'calib_IO_Cam2.mat']; %thermal to RGB transformation matrix
fileCalibStereo_RGB2TIR = [WORKING_FOLDER 'myRbTb_fromCalTech3D.mat']; %results from CalTech Camera Calib Stereo RGB to theraml bore sight

% mask file
MASK_FILE = [WORKING_FOLDER 'mask.png'];

%% %%%%%%%%%%%%%%%%%%%%%%%%% GET CAMERA PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ CAMERA EO/IO AND CALIBRATION
cameraXML = xml2struct(fileCameraRGB_XML); % read camera EO and IO from an XML file to a structure
camTermalCalib = load(fileThermalCameraCalib, '-mat'); % read thermal camera params
boreSightRGB2TIR = load(fileCalibStereo_RGB2TIR, '-mat'); % read RGB to TIR transformation

% GET PARAMS
camIO = GetCamIO(cameraXML);  % get Interior Orinetation Parameters
camEO = GetCamEO(cameraXML); % get Exterior Orinetation Parameters
[Rotation, Translation, Scale] = GetWorldTransform(cameraXML);  % get Coordiante Transforamtion Matrix

nCam = size(cameraXML.document.chunk.cameras.camera,2); % number of images

% read mask png file
MASK = imread(MASK_FILE);

%% %%%%%%%%%%%%%%%%%%%%%%%% READ TIR IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiate read image waitbar
readImgWaintbar = waitbar(0,'Reading TIR images...');


% read in thermal images and store them in a cell array
imArrayTIR = cell(1,nCam); %preallocate a struct to store TIR images
for iCam=1:nCam
    imageFileName = [FOLDER_TIR_IMAGES camEO(iCam).label]; %get the RGB image name (.jpg)
    imageNameTIR = [imageFileName(1:end-4) '.jpg']; %generate thermal image name (.tif)
    
    img_Thermal = imread(imageNameTIR);
    
    if (TO_APPLY_MASK == 1)
        img_Thermal = AddMaskToImage(img_Thermal, MASK);
    end
    
    imArrayTIR{iCam} = img_Thermal(:,:,1); %stote all images in an image array
    
    % update read image waitbar
    waitbar(iCam/nCam,readImgWaintbar,'Reading TIR images...','CreateCancelBtn','button_callback');
end

%get min and max values on thermal images
[minT, maxT] = GetMinMaxFromStruct(imArrayTIR);

% close read image waitbar
delete(readImgWaintbar);


%% %%%%%%%%%%%%%%%%%%%%%%%%% READ POINTCLOUD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiate read file waitbar
readCloudWaitbar = waitbar(0,'Reading point cloud...');

% Read point cloud
cloudData = load(fileCloudRead);
% number if points
nPts = size(cloudData, 1);

% clode read file waitbar
delete(readCloudWaitbar);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open a pointcloud file and write
fileID = fopen(fileCloudWrite,'w');

% Update waitbar
generateCloudWaitbar =waitbar(0,'Generating thermal point cloud...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(generateCloudWaitbar,'canceling',0)

for pt = 1:nPts
    % Check for Cancel button press
    if getappdata(generateCloudWaitbar,'canceling') end
    
    % preallocate a variable to store thermal reading from images
    readingsThermalValue = NaN(1,nCam);
    
    %read world coordiante of the point
    
    XYZ_w = transpose(cloudData(pt,1:3));
    
    %get chuck coordinate
    xyz_block = TransWorldToBlock(XYZ_w, Rotation, Translation, Scale);
    
    for iCam = 1: nCam
        
        %get coordinate of the point with respect to camera iCam
        xyz_c = TransBlockToCam(xyz_block, camEO(iCam));
        
        %apply boresight traansforamtin to get the coordiantes the point with respect to TIR camera using
        xyz_cT = TransBoresight(xyz_c, boreSightRGB2TIR);
        
        % find matching pixel in thermal image
        % uT (right) vT (down)
        [uT, vT] = TransCamToImage_TIR(xyz_cT, camTermalCalib.camCalibIO);
        
        if (IsInImage(uT,vT, 2*camTermalCalib.camCalibIO.cc(2), 2*camTermalCalib.camCalibIO.cc(1) ) == true) %in thermal
            
            % get the thermal value
            img_TIR = imArrayTIR{iCam};
            readingsThermalValue(iCam) = GetPixelValue(img_TIR, uT, vT, 1);
            
            %             if TO_DISPLAY == 1
            %                 % calculate image coordiantes of the point
            %                 [u, v] = TransCamToImage(xyz_c, camIO);
            %                 img_RGB = imread(strcat(FOLDER_RGB_IMAGES, camEO(iCam).label));
            %                 DisplayAndPlot(img_RGB, img_TIR, u, v, uT, vT);
            %             end %if display
            
        end %if in image
        
    end %for iCam
    intensity = round(nanmedian(readingsThermalValue)); %get integer value
    
    %convert intensity to absolute thermal value
    tempC = ConvertIntesityToTempFrom_MIN_MAX(intensity, minT, maxT, TEMP_MIN_READING,TEMP_MAX_READING );
    
    
    %write out to file line by line
    fprintf(fileID, '%.6f %.6f %.6f %d %d %d %0.2f\n',...
        XYZ_w(1), XYZ_w(2), XYZ_w(3), cloudData(pt,4), cloudData(pt,5), cloudData(pt,6), tempC);
    
    % Update waitbar
    if round(pt/nPts*100) > round((pt-1)/nPts*100)
        waitbar(pt/nPts,generateCloudWaitbar, sprintf('Generating thermal point cloud... (%0.0f%%)',100*pt/nPts));
    end
end %for pt

% close the file
fclose(fileID);

%clode waitbar
delete(generateCloudWaitbar);
