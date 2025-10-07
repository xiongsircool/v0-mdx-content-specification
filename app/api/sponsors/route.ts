import { NextResponse } from 'next/server'
import fs from 'fs'
import path from 'path'

export async function GET() {
  try {
    const sponsorsPath = path.join(process.cwd(), 'content', 'sponsors', 'sponsors.json')
    const sponsorsData = fs.readFileSync(sponsorsPath, 'utf-8')
    const sponsors = JSON.parse(sponsorsData)

    return NextResponse.json(sponsors)
  } catch (error) {
    console.error('Error loading sponsors:', error)
    return NextResponse.json([], { status: 500 })
  }
}