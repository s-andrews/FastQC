<?xml version="1.0"?>
<xsl:stylesheet 
  version="1.0" 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:fo="http://www.w3.org/1999/XSL/Format"
  xmlns:fox="http://xml.apache.org/fop/extensions">

  <xsl:template match="html">
    <fo:root xmlns:fo="http://www.w3.org/1999/XSL/Format"
      xmlns:fox="http://xml.apache.org/fop/extensions"> 

     <fo:layout-master-set>
        <fo:simple-page-master master-name="page-layout">
          <fo:region-body margin="2.5cm" region-name="body"/>
        </fo:simple-page-master>
      </fo:layout-master-set>

      <fo:page-sequence master-reference="page-layout">



    <xsl:apply-templates select="body"/>
      </fo:page-sequence>
    </fo:root>
  </xsl:template>


  <xsl:template match="b">
    <fo:inline font-weight="bold">
      <xsl:apply-templates select="*|text()"/>
    </fo:inline>
  </xsl:template>




  <xsl:template match="body">
  
    <fo:flow flow-name="body">
        <fo:block font-size="48pt" text-align="center">
              FASTQC-Report
            <fo:inline wrap-option="no-wrap"/>
         </fo:block>
 		 
           <xsl:apply-templates select="//div[@class='module']"/>

         </fo:flow>
  </xsl:template>

<xsl:template match="div[@class='module']">
     <fo:block text-align="center" page-break-before="always"  font-size="48pt"> 
           <xsl:value-of select="h2"/>
      </fo:block>
	<xsl:apply-templates select="p/img|table"/>

</xsl:template>

  <xsl:template match="img">

    <xsl:if  test="starts-with(@src,'Images/')">
        <fo:block space-after="12pt" width="18cm" >
          <fo:external-graphic src="{@src}"  content-width="19cm" content-height="12.4cm" scaling="uniform" >
          </fo:external-graphic>
        </fo:block>
    </xsl:if>
  </xsl:template>



  <xsl:template match="table">
    <fo:table >
        <xsl:for-each select="thead/tr/th">
      		<fo:table-column column-width="200pt"/>
        </xsl:for-each>
      <fo:table-body>
        <xsl:apply-templates select="tbody/tr"/>
      </fo:table-body>
    </fo:table>
  </xsl:template>

  <xsl:template match="tr">
    <fo:table-row>
      <xsl:apply-templates select="td"/>
    </fo:table-row>
  </xsl:template>

  <xsl:template match="td">
    <fo:table-cell>
      <fo:block>
        <xsl:apply-templates select="*|text()"/>
      </fo:block>
    </fo:table-cell>
  </xsl:template>


</xsl:stylesheet>