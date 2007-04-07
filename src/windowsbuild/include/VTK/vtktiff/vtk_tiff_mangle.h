#ifndef vtk_tiff_mangle_h
#define vtk_tiff_mangle_h

/*

This header file mangles all symbols exported from the tiff library.
It is included in all files while building the tiff library.  Due to
namespace pollution, no tiff headers should be included in .h files in
VTK.

The following command was used to obtain the symbol list:

nm libvtktiff.a |grep " [TR] "

*/

#define __TIFFmalloc vtk___TIFFmalloc
#define __TIFFfree vtk___TIFFfree
#define _TIFFerrorHandler vtk__TIFFerrorHandler
#define _TIFFwarningHandler vtk__TIFFwarningHandler
#define tiffDataWidth vtk_tiff_DataWidth
#define LogL10fromY vtk_tiff_LogL10fromY
#define LogL10toY vtk_tiff_LogL10toY
#define LogL16fromY vtk_tiff_LogL16fromY
#define LogL16toY vtk_tiff_LogL16toY
#define LogLuv24fromXYZ vtk_tiff_LogLuv24fromXYZ
#define LogLuv24toXYZ vtk_tiff_LogLuv24toXYZ
#define LogLuv32fromXYZ vtk_tiff_LogLuv32fromXYZ
#define LogLuv32toXYZ vtk_tiff_LogLuv32toXYZ
#define TIFFCheckTile vtk_TIFFCheckTile
#define TIFFClientOpen vtk_TIFFClientOpen
#define TIFFClose vtk_TIFFClose
#define TIFFComputeStrip vtk_TIFFComputeStrip
#define TIFFComputeTile vtk_TIFFComputeTile
#define TIFFCreateDirectory vtk_TIFFCreateDirectory
#define TIFFCurrentDirOffset vtk_TIFFCurrentDirOffset
#define TIFFCurrentDirectory vtk_TIFFCurrentDirectory
#define TIFFCurrentRow vtk_TIFFCurrentRow
#define TIFFCurrentStrip vtk_TIFFCurrentStrip
#define TIFFCurrentTile vtk_TIFFCurrentTile
#define TIFFDefaultDirectory vtk_TIFFDefaultDirectory
#define TIFFDefaultStripSize vtk_TIFFDefaultStripSize
#define TIFFDefaultTileSize vtk_TIFFDefaultTileSize
#define TIFFError vtk_TIFFError
#define TIFFFaxBlackCodes vtk_TIFFFaxBlackCodes
#define TIFFFaxBlackTable vtk_TIFFFaxBlackTable
#define TIFFFaxMainTable vtk_TIFFFaxMainTable
#define TIFFFaxWhiteCodes vtk_TIFFFaxWhiteCodes
#define TIFFFaxWhiteTable vtk_TIFFFaxWhiteTable
#define TIFFFdOpen vtk_TIFFFdOpen
#define TIFFFileName vtk_TIFFFileName
#define TIFFFileno vtk_TIFFFileno
#define TIFFFindCODEC vtk_TIFFFindCODEC
#define TIFFFlush vtk_TIFFFlush
#define TIFFFlushData vtk_TIFFFlushData
#define TIFFFlushData1 vtk_TIFFFlushData1
#define TIFFFreeDirectory vtk_TIFFFreeDirectory
#define TIFFGetBitRevTable vtk_TIFFGetBitRevTable
#define TIFFGetField vtk_TIFFGetField
#define TIFFGetFieldDefaulted vtk_TIFFGetFieldDefaulted
#define TIFFGetMode vtk_TIFFGetMode
#define TIFFGetVersion vtk_TIFFGetVersion
#define TIFFInitCCITTFax3 vtk_TIFFInitCCITTFax3
#define TIFFInitCCITTFax4 vtk_TIFFInitCCITTFax4
#define TIFFInitCCITTRLE vtk_TIFFInitCCITTRLE
#define TIFFInitCCITTRLEW vtk_TIFFInitCCITTRLEW
#define TIFFInitDumpMode vtk_TIFFInitDumpMode
#define TIFFInitJPEG vtk_TIFFInitJPEG
#define TIFFInitLZW vtk_TIFFInitLZW
#define TIFFInitNeXT vtk_TIFFInitNeXT
#define TIFFInitPackBits vtk_TIFFInitPackBits
#define TIFFInitSGILog vtk_TIFFInitSGILog
#define TIFFInitThunderScan vtk_TIFFInitThunderScan
#define TIFFInitZIP vtk_TIFFInitZIP
#define TIFFIsByteSwapped vtk_TIFFIsByteSwapped
#define TIFFIsMSB2LSB vtk_TIFFIsMSB2LSB
#define TIFFIsTiled vtk_TIFFIsTiled
#define TIFFIsUpSampled vtk_TIFFIsUpSampled
#define TIFFLastDirectory vtk_TIFFLastDirectory
#define TIFFNumberOfDirectories vtk_TIFFNumberOfDirectories
#define TIFFNumberOfStrips vtk_TIFFNumberOfStrips
#define TIFFNumberOfTiles vtk_TIFFNumberOfTiles
#define TIFFOpen vtk_TIFFOpen
#define TIFFPredictorInit vtk_TIFFPredictorInit
#define TIFFPrintDirectory vtk_TIFFPrintDirectory
#define TIFFRGBAImageBegin vtk_TIFFRGBAImageBegin
#define TIFFRGBAImageEnd vtk_TIFFRGBAImageEnd
#define TIFFRGBAImageGet vtk_TIFFRGBAImageGet
#define TIFFRGBAImageOK vtk_TIFFRGBAImageOK
#define TIFFRasterScanlineSize vtk_TIFFRasterScanlineSize
#define TIFFReadBufferSetup vtk_TIFFReadBufferSetup
#define TIFFReadDirectory vtk_TIFFReadDirectory
#define TIFFReadEncodedStrip vtk_TIFFReadEncodedStrip
#define TIFFReadEncodedTile vtk_TIFFReadEncodedTile
#define TIFFReadRGBAImage vtk_TIFFReadRGBAImage
#define TIFFReadRGBAStrip vtk_TIFFReadRGBAStrip
#define TIFFReadRGBATile vtk_TIFFReadRGBATile
#define TIFFReadRawStrip vtk_TIFFReadRawStrip
#define TIFFReadRawTile vtk_TIFFReadRawTile
#define TIFFReadScanline vtk_TIFFReadScanline
#define TIFFReadTile vtk_TIFFReadTile
#define TIFFReassignTagToIgnore vtk_TIFFReassignTagToIgnore
#define TIFFRegisterCODEC vtk_TIFFRegisterCODEC
#define TIFFReverseBits vtk_TIFFReverseBits
#define TIFFRewriteDirectory vtk_TIFFRewriteDirectory
#define TIFFScanlineSize vtk_TIFFScanlineSize
#define TIFFSetCompressionScheme vtk_TIFFSetCompressionScheme
#define TIFFSetDirectory vtk_TIFFSetDirectory
#define TIFFSetErrorHandler vtk_TIFFSetErrorHandler
#define TIFFSetField vtk_TIFFSetField
#define TIFFSetSubDirectory vtk_TIFFSetSubDirectory
#define TIFFSetTagExtender vtk_TIFFSetTagExtender
#define TIFFSetWarningHandler vtk_TIFFSetWarningHandler
#define TIFFSetWriteOffset vtk_TIFFSetWriteOffset
#define TIFFStripSize vtk_TIFFStripSize
#define TIFFSwabArrayOfDouble vtk_TIFFSwabArrayOfDouble
#define TIFFSwabArrayOfLong vtk_TIFFSwabArrayOfLong
#define TIFFSwabArrayOfShort vtk_TIFFSwabArrayOfShort
#define TIFFSwabDouble vtk_TIFFSwabDouble
#define TIFFSwabLong vtk_TIFFSwabLong
#define TIFFSwabShort vtk_TIFFSwabShort
#define TIFFTileRowSize vtk_TIFFTileRowSize
#define TIFFTileSize vtk_TIFFTileSize
#define TIFFUnRegisterCODEC vtk_TIFFUnRegisterCODEC
#define TIFFUnlinkDirectory vtk_TIFFUnlinkDirectory
#define TIFFVGetField vtk_TIFFVGetField
#define TIFFVGetFieldDefaulted vtk_TIFFVGetFieldDefaulted
#define TIFFVSetField vtk_TIFFVSetField
#define TIFFVStripSize vtk_TIFFVStripSize
#define TIFFVTileSize vtk_TIFFVTileSize
#define TIFFWarning vtk_TIFFWarning
#define TIFFWriteBufferSetup vtk_TIFFWriteBufferSetup
#define TIFFWriteCheck vtk_TIFFWriteCheck
#define TIFFWriteDirectory vtk_TIFFWriteDirectory
#define TIFFWriteEncodedStrip vtk_TIFFWriteEncodedStrip
#define TIFFWriteEncodedTile vtk_TIFFWriteEncodedTile
#define TIFFWriteRawStrip vtk_TIFFWriteRawStrip
#define TIFFWriteRawTile vtk_TIFFWriteRawTile
#define TIFFWriteScanline vtk_TIFFWriteScanline
#define TIFFWriteTile vtk_TIFFWriteTile
#define XYZtoRGB24 vtk_tiff_XYZtoRGB24
#define _TIFFDefaultStripSize vtk__TIFFDefaultStripSize
#define _TIFFDefaultTileSize vtk__TIFFDefaultTileSize
#define _TIFFFax3fillruns vtk__TIFFFax3fillruns
#define _TIFFFieldWithTag vtk__TIFFFieldWithTag
#define _TIFFFindFieldInfo vtk__TIFFFindFieldInfo
#define _TIFFMergeFieldInfo vtk__TIFFMergeFieldInfo
#define _TIFFNoPostDecode vtk__TIFFNoPostDecode
#define _TIFFNoPreCode vtk__TIFFNoPreCode
#define _TIFFNoRowDecode vtk__TIFFNoRowDecode
#define _TIFFNoRowEncode vtk__TIFFNoRowEncode
#define _TIFFNoSeek vtk__TIFFNoSeek
#define _TIFFNoStripDecode vtk__TIFFNoStripDecode
#define _TIFFNoStripEncode vtk__TIFFNoStripEncode
#define _TIFFNoTileDecode vtk__TIFFNoTileDecode
#define _TIFFNoTileEncode vtk__TIFFNoTileEncode
#define _TIFFPrintFieldInfo vtk__TIFFPrintFieldInfo
#define _TIFFSampleToTagType vtk__TIFFSampleToTagType
#define _TIFFSetDefaultCompressionState vtk__TIFFSetDefaultCompressionState
#define _TIFFSetupFieldInfo vtk__TIFFSetupFieldInfo
#define _TIFFSwab16BitData vtk__TIFFSwab16BitData
#define _TIFFSwab32BitData vtk__TIFFSwab32BitData
#define _TIFFSwab64BitData vtk__TIFFSwab64BitData
#define _TIFFfree vtk__TIFFfree
#define _TIFFgetMode vtk__TIFFgetMode
#define _TIFFmalloc vtk__TIFFmalloc
#define _TIFFmemcmp vtk__TIFFmemcmp
#define _TIFFmemcpy vtk__TIFFmemcpy
#define _TIFFmemset vtk__TIFFmemset
#define _TIFFprintAscii vtk__TIFFprintAscii
#define _TIFFprintAsciiTag vtk__TIFFprintAsciiTag
#define _TIFFrealloc vtk__TIFFrealloc
#define _TIFFsetByteArray vtk__TIFFsetByteArray
#define _TIFFsetDoubleArray vtk__TIFFsetDoubleArray
#define _TIFFsetFloatArray vtk__TIFFsetFloatArray
#define _TIFFsetLongArray vtk__TIFFsetLongArray
#define _TIFFsetNString vtk__TIFFsetNString
#define _TIFFsetShortArray vtk__TIFFsetShortArray
#define _TIFFsetString vtk__TIFFsetString
#define uv_decode vtk_tiff_uv_decode
#define uv_encode vtk_tiff_uv_encode

#endif
