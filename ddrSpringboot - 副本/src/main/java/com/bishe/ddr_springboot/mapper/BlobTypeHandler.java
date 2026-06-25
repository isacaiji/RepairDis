package com.bishe.ddr_springboot.mapper;

import org.apache.ibatis.type.BaseTypeHandler;
import org.apache.ibatis.type.JdbcType;
import org.apache.ibatis.type.MappedJdbcTypes;
import org.apache.ibatis.type.MappedTypes;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.sql.CallableStatement;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

@MappedJdbcTypes(JdbcType.BLOB)
@MappedTypes(String.class)
public class BlobTypeHandler extends BaseTypeHandler<String> {
    private static final Logger logger = LoggerFactory.getLogger(BlobTypeHandler.class);

    @Override
    public void setNonNullParameter(PreparedStatement ps, int i, String parameter, JdbcType jdbcType) throws SQLException {
        try {
            InputStream inputStream = new ByteArrayInputStream(parameter.getBytes("UTF-8"));
            ps.setBlob(i, inputStream);
        } catch (UnsupportedEncodingException e) {
            logger.error("Error setting BLOB parameter", e);
            throw new RuntimeException(e);
        }
    }

    @Override
    public String getNullableResult(ResultSet rs, String columnName) throws SQLException {
        java.sql.Blob blob = rs.getBlob(columnName);
        if (blob != null) {
            logger.info("Trying to convert BLOB for column: {}", columnName);
            try (InputStream inputStream = blob.getBinaryStream()) {
                byte[] bytes = new byte[(int) blob.length()];
                int readBytes = inputStream.read(bytes);
                logger.info("Read {} bytes from BLOB for column: {}", readBytes, columnName);
                return new String(bytes, "UTF-8");
            } catch (IOException e) {
                logger.error("Error reading BLOB data for column: {}", columnName, e);
                return null;
            }
        } else {
            logger.info("BLOB value is null for column: {}", columnName);
        }
        return null;
    }

    @Override
    public String getNullableResult(ResultSet rs, int columnIndex) throws SQLException {
        java.sql.Blob blob = rs.getBlob(columnIndex);
        if (blob != null) {
            logger.info("Trying to convert BLOB for column index: {}", columnIndex);
            try (InputStream inputStream = blob.getBinaryStream()) {
                byte[] bytes = new byte[(int) blob.length()];
                int readBytes = inputStream.read(bytes);
                logger.info("Read {} bytes from BLOB for column index: {}", readBytes, columnIndex);
                return new String(bytes, "UTF-8");
            } catch (IOException e) {
                logger.error("Error reading BLOB data for column index: {}", columnIndex, e);
                return null;
            }
        } else {
            logger.info("BLOB value is null for column index: {}", columnIndex);
        }
        return null;
    }
    @Override
    public String getNullableResult(CallableStatement cs, int columnIndex) throws SQLException {
        java.sql.Blob blob = cs.getBlob(columnIndex);
        if (blob != null) {
            logger.info("Trying to convert BLOB for column index: {}", columnIndex);
            try (InputStream inputStream = blob.getBinaryStream()) {
                byte[] bytes = new byte[(int) blob.length()];
                int readBytes = inputStream.read(bytes);
                logger.info("Read {} bytes from BLOB for column index: {}", readBytes, columnIndex);
                return new String(bytes, "UTF-8");
            } catch (IOException e) {
                logger.error("Error reading BLOB data for column index: {}", columnIndex, e);
                return null;
            }
        } else {
            logger.info("BLOB value is null for column index: {}", columnIndex);
        }
        return null;
    }
}