#include "mex.h"
#include <stdio.h>
#include <Windows.h>
#include <MSR_NuiApi.h>
#define imageHeight 480
#define imageWidth  640
bool Nui_GotDepthAlert(double  Depthes[][imageWidth]);
bool Nui_GotVideoAlert(double  RGBinfo[][imageWidth][3]);
void KinectInit();
void KinectShutdown();
HRESULT	hr;
HANDLE m_pDepthStreamHandle = NULL;
HANDLE m_pVideoStreamHandle;
const NUI_IMAGE_FRAME * pImageFrame = NULL;

int main(){}

void KinectInit()
{
	hr=NuiInitialize( NUI_INITIALIZE_FLAG_USES_DEPTH | NUI_INITIALIZE_FLAG_USES_COLOR);

	hr = NuiImageStreamOpen(NUI_IMAGE_TYPE_COLOR,
         NUI_IMAGE_RESOLUTION_640x480, 0, 2, NULL, &m_pVideoStreamHandle );

	hr = NuiImageStreamOpen( NUI_IMAGE_TYPE_DEPTH,
		 NUI_IMAGE_RESOLUTION_640x480, 0, 2, NULL, &m_pDepthStreamHandle);
}

void KinectShutdown()
{
	NuiShutdown();
}

bool Nui_GotDepthAlert(double Depthes[][imageWidth])
{
    const NUI_IMAGE_FRAME * pImageFrame = NULL;
    HRESULT hr = NuiImageStreamGetNextFrame(
        m_pDepthStreamHandle,	
        66,//2/30 millisecond timeout
        &pImageFrame );
    if( FAILED( hr ) )
		return -1;
    NuiImageBuffer * pTexture = pImageFrame->pFrameTexture;
    KINECT_LOCKED_RECT LockedRect;
    pTexture->LockRect( 0, &LockedRect, NULL, 0 );
    if( LockedRect.Pitch != 0 )
    {
        BYTE * pBuffer = (BYTE*) LockedRect.pBits;
        USHORT * pBufferRun = (USHORT*) pBuffer;
		double RealDepth;
		for( int y = 0 ; y < imageHeight ; y++ )
        {	
            for( int x = 0 ; x < imageWidth ; x++ )
            {
                pBufferRun++;
				RealDepth = (*pBufferRun  & 0x0fff);
				Depthes[y][x] = RealDepth;		
            }
        }	
    }
    else
    {
	  //OutputDebugString( L"Buffer length of received texture is bogus\r\n" );
    }
    NuiImageStreamReleaseFrame( m_pDepthStreamHandle, pImageFrame );
	return 0;
}

bool Nui_GotVideoAlert(double RGBinfo[][imageWidth][3])
{
    const NUI_IMAGE_FRAME * pImageFrame = NULL;

    HRESULT hr = NuiImageStreamGetNextFrame(
        m_pVideoStreamHandle,
        66,
        &pImageFrame );
    if( FAILED( hr ) )
        return -1;
    NuiImageBuffer * pTexture = pImageFrame->pFrameTexture;
    KINECT_LOCKED_RECT LockedRect;
    pTexture->LockRect( 0, &LockedRect, NULL, 0 );
    if( LockedRect.Pitch != 0 )
    {
        BYTE * pBuffer = (BYTE*) LockedRect.pBits;
		int base, b, a= 0;
		for (int j= 0;j< imageHeight*4;j+= 4){
			b= 0;
			base = j*imageWidth;
			for(int i= 0;i< (imageWidth*4);i+= 4)
			{
				RGBinfo[a][b][0]= (double)pBuffer[base+i+2]; //R
				RGBinfo[a][b][1]= (double)pBuffer[base+i+1]; //G
				RGBinfo[a][b][2]= (double)pBuffer[base+i+0]; //B
				b++;
			}
			a++;
		}
    }
	else
    {
		//OutputDebugString( L"Buffer length of received texture is bogus\r\n" ); 
    }
    NuiImageStreamReleaseFrame( m_pVideoStreamHandle, pImageFrame );
	return 0;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{	
	    if(nlhs!=2)
		{
			mexErrMsgIdAndTxt("Getimagedata:nlhs","One output required.");
		}
		double *Dep, *Rgb, *Option;
		int dims[3] = {480, 640, 3};
		plhs[0] = mxCreateDoubleMatrix(480, 640, mxREAL);
		plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		Dep = mxGetPr(plhs[0]);
		Rgb = mxGetPr(plhs[1]);
		Option = mxGetPr(prhs[0]);
		double (*Depxy)[640] = (double (*)[640])Dep;
		double (*Rgbxyz)[640][3] = (double (*)[640][3])Rgb;
		switch ((int)*Option)
		{
			case 1 :
				KinectInit();			
				mexLock();
				break;
			case 2 :
				Nui_GotDepthAlert(Depxy);
				Nui_GotVideoAlert(Rgbxyz);;
				break;
			case 3:
				/*
				KinectShutdown();
				mexUnlock();
				*/
				break;
			default:
				break;
		}
}